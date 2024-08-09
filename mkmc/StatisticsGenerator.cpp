#include "StatisticsGenerator.h"
#include "MatrixStats.h"
#include "DumpWriter.h"
#include <algorithm>



void StatisticsGenerator::openReaders()
{
	try
	{
		matrixMetadataReader = std::make_unique<kmcdb::MetadataReader>(params.mkmcParams.outputBinFile, false);
		matrixReader = std::make_unique<kmcdb::ReaderSortedPlainForListing<uint64_t>>(*matrixMetadataReader);
	}
	catch (const std::runtime_error& ex)
	{
		std::cerr << "Error: " << ex.what() << std::endl;
		exit(1);
	}
}



void StatisticsGenerator::fillTaskData()
{
	tasksData.reserve(params.stage1Params.GetNBins());
	if (params.statisticsParams.correctPvalues)
		correctTasksData.reserve(statisticsToGeneration.nStatisticsWithPValues);

	nOutputKmersPerBin.reserve(params.stage1Params.GetNBins());
	
	for (uint32_t i = 0; i < params.stage1Params.GetNBins(); ++i)
	{
		tasksData.push_back(TaskData{ i });
		nOutputKmersPerBin.push_back(matrixReader->GetBin(i)->GetBinMetadata().total_kmers);
	}

	if (params.statisticsParams.correctPvalues)
		for (uint32_t i = 0; i < statisticsToGeneration.nStatisticsWithPValues; ++i)
			correctTasksData.push_back(CorrectTaskData{ i });

	for (size_t it = 1; it < binsIndicesForCorrection.size(); ++it) // fragmentsBegins[0] = 0
	{
		binsIndicesForCorrection[it] = binsIndicesForCorrection[it - 1] + nOutputKmersPerBin[it - 1]; // increase previous index by a size of the next bin
	}

	uint64_t progressBarTicks = params.mkmcParams.verbosity_level == 0 ? 0 : std::accumulate(nOutputKmersPerBin.begin(), nOutputKmersPerBin.end(), 0ull);
	if (params.statisticsParams.correctPvalues)
		progressBarTicks *= 2;
	progress_bar = std::make_unique<ProgressBar>(
		progressBarTicks,
		"Computing statistics",
		std::cerr,
		params.mkmcParams.verbosity_level == 0);

	std::sort(tasksData.begin(), tasksData.end(), [&](const TaskData& a, const TaskData& b) { return nOutputKmersPerBin[a.binId] > nOutputKmersPerBin[b.binId]; });
	// correctTasksData does need to be sorted
}



StatisticsGenerator::StatisticsGenerator(Params& params) :
	params(params),
	tasksPool(tasksData),
	correctTasksPool(correctTasksData),
	gatherer(params, statisticsToGeneration),
	correlationPhenotype(params.phenotypes.correlationPhenotype.getPhenotype()),
	differentialAnalysisPhenotype(params.phenotypes.differentialAnalysisPhenotype.getMappedPhenotype()),
	differentialAnalysisNClasses(params.phenotypes.differentialAnalysisPhenotype.getClassesNumber()),
	binsIndicesForCorrection(params.stage1Params.GetNBins() + 1)
{

	auto is_correlation_method = [&](StatisticsParams::CorrelationMethod method)
	{
		const auto& corMeths = params.statisticsParams.correlationMethods;
		if (std::find(corMeths.begin(), corMeths.end(), method) != corMeths.end())
		{
			++statisticsToGeneration.nStatistics;
			return true;
		}
		return false;
	};

	auto is_differential_analysis_method = [&](StatisticsParams::DifferentialAnalysisMethod method)
	{
		const auto& analysisMeths = params.statisticsParams.classificationMethods;
		if (std::find(analysisMeths.begin(), analysisMeths.end(), method) != analysisMeths.end())
		{
			++statisticsToGeneration.nStatistics;
			return true;
		}
		return false;
	};

	statisticsToGeneration.normalize = params.statisticsParams.generateNormalization;

	statisticsToGeneration.pearson = is_correlation_method(StatisticsParams::CorrelationMethod::Pearson);
	statisticsToGeneration.spearman = is_correlation_method(StatisticsParams::CorrelationMethod::Spearman);
	statisticsToGeneration.kendall = is_correlation_method(StatisticsParams::CorrelationMethod::Kendall);

	statisticsToGeneration.entropy = params.statisticsParams.generateEntropy;
	if (statisticsToGeneration.entropy)
		++statisticsToGeneration.nStatistics;

	statisticsToGeneration.tTest = is_differential_analysis_method(StatisticsParams::DifferentialAnalysisMethod::TTest);
	statisticsToGeneration.snr = is_differential_analysis_method(StatisticsParams::DifferentialAnalysisMethod::SNR);
	statisticsToGeneration.wilcoxonRankSum = is_differential_analysis_method(StatisticsParams::DifferentialAnalysisMethod::WilcoxonRankSum);

	statisticsToGeneration.dids = is_differential_analysis_method(StatisticsParams::DifferentialAnalysisMethod::DIDS);
	statisticsToGeneration.anova = is_differential_analysis_method(StatisticsParams::DifferentialAnalysisMethod::ANOVA);

	statisticsToGeneration.differentialAnalysis = statisticsToGeneration.tTest || statisticsToGeneration.snr || statisticsToGeneration.wilcoxonRankSum || statisticsToGeneration.dids || statisticsToGeneration.anova;

	statisticsToGeneration.nResults = statisticsToGeneration.nStatistics + (params.statisticsParams.generateNormalization ? params.mkmcParams.samples.size() : 0);

	if (statisticsToGeneration.tTest)
		++statisticsToGeneration.nStatisticsWithPValues;
	if (statisticsToGeneration.wilcoxonRankSum)
		++statisticsToGeneration.nStatisticsWithPValues;
	if (statisticsToGeneration.anova)
		++statisticsToGeneration.nStatisticsWithPValues;
}



void StatisticsGenerator::processEntries()
{
	TaskData taskData;
	while (tasksPool.getTask(taskData))
	{
		auto bin = matrixReader->GetBin(taskData.binId);
		std::unique_ptr<OutputBuffer> normOutputBuffer;

		refresh::normalization_work<uint64_t, double> normalization;
		if (params.statisticsParams.generateNormalization)
		{
			normalization.register_method(params.statisticsParams.normalizationMethod);
			normalization.set_no_series(params.mkmcParams.samples.size());
			normalization.deserialize(params.statisticsParams.normalizationMethod, normalizationData);

			normalization.initialize();

			normOutputBuffer = std::make_unique<OutputBuffer>(*normWriter, getMaxNormLineLength());
		}

		refresh::correlation correlation;
		refresh::statistics_entropy entropyObj;
		refresh::statistical_test statistics;
		refresh::scorers scorer;

		std::vector<uint64_t> inMatrixEntry;
		std::vector<out_kmcdb_value_type> outEntry; // normalized counts and statistics
		std::ptrdiff_t num_samples = static_cast<std::ptrdiff_t>(params.mkmcParams.samples.size());

		inMatrixEntry.resize(num_samples);
		outEntry.resize(statisticsToGeneration.nResults);

		ProgressBarUpdater progress_bar_updater(*progress_bar, (std::max)(1ull, progress_bar->GetTotal() / 100ull));

		auto kmer_len = params.stage1Params.GetKmerLen();
		std::string kmerSequence(kmer_len, ' ');

		kmcdb::DispatchKmerSize<MAX_K>(kmer_len, [&](auto SIZE) {
			kmcdb::CKmer<SIZE> kmer;
			std::unique_ptr<WritingGathererBin<out_kmcdb_value_type>> outBin = gatherer.getBin(taskData.binId);
			while (bin->NextKmer(kmer, inMatrixEntry.data()))
			{
				kmer.to_string(kmer_len, kmerSequence.data());

				if (statisticsToGeneration.normalize)
				{
					normalization.norm_entry(params.statisticsParams.normalizationMethod, inMatrixEntry, outEntry);

					normOutputBuffer->StoreKmer(kmerSequence, outEntry, StoreMethods::AsMatrixRow);

					outEntry.resize(statisticsToGeneration.nResults); // space for statistics
				}

				size_t outStatsEntryIdx = statisticsToGeneration.nResults - statisticsToGeneration.nStatistics;

				if (statisticsToGeneration.pearson)
				{
					assert(statisticsToGeneration.normalize);
					const double pearson = refresh::correlation::pearson_n(
						outEntry.begin(),
						correlationPhenotype.begin(),
						num_samples);

					outEntry[outStatsEntryIdx++] = pearson;
				}
				if (statisticsToGeneration.spearman)
				{
					assert(statisticsToGeneration.normalize);
					const double spearman = correlation.spearman_n(
						outEntry.begin(),
						correlationPhenotype.begin(),
						num_samples);

					outEntry[outStatsEntryIdx++] = spearman;
				}
				if (statisticsToGeneration.kendall)
				{
					assert(statisticsToGeneration.normalize);
					const double kendall = refresh::correlation::kendall_tau_n(
						outEntry.begin(),
						correlationPhenotype.begin(),
						num_samples);

					outEntry[outStatsEntryIdx++] = kendall;
				}

				if (statisticsToGeneration.entropy)
				{
					const double entropy = entropyObj.entropy_n(
						inMatrixEntry.begin(),
						num_samples);

					outEntry[outStatsEntryIdx++] = entropy;
				}
				if (statisticsToGeneration.differentialAnalysis)
				{
					if (statisticsToGeneration.tTest)
					{
						const double tTestPValue = statistics.t_test_n(
							inMatrixEntry.begin(),
							differentialAnalysisPhenotype.begin(),
							num_samples).p_value;

						outEntry[outStatsEntryIdx++] = tTestPValue;
					}
					if (statisticsToGeneration.snr)
					{
						const double snr = statistics.SNR_test_n(
							inMatrixEntry.begin(),
							differentialAnalysisPhenotype.begin(),
							num_samples);

						outEntry[outStatsEntryIdx++] = snr;
					}
					if (statisticsToGeneration.wilcoxonRankSum)
					{
						const double wilcoxonRankSumPValue = statistics.mann_whitney_U_test_n(
							inMatrixEntry.begin(),
							differentialAnalysisPhenotype.begin(),
							num_samples).p_value;

						outEntry[outStatsEntryIdx++] = wilcoxonRankSumPValue;
					}
					if (statisticsToGeneration.dids)
					{
						const double dids = scorer.dids_n(
							inMatrixEntry.begin(),
							differentialAnalysisPhenotype.begin(),
							differentialAnalysisNClasses,
							num_samples);

						outEntry[outStatsEntryIdx++] = dids;
					}
					if (statisticsToGeneration.anova)
					{
						const double anovaPValue = scorer.anova_n(
							inMatrixEntry.begin(),
							differentialAnalysisPhenotype.begin(),
							differentialAnalysisNClasses,
							num_samples).p_value;

						outEntry[outStatsEntryIdx++] = anovaPValue;
					}
				}

				++progress_bar_updater;

				outBin->writeKmer(outEntry, kmer, kmerSequence, inMatrixEntry);
			}
		});
	}
}



void StatisticsGenerator::gatherPValuesEntriesToCorrection()
{
	assert(statisticsToGeneration.differentialAnalysis);
	TaskData taskData;
	while (tasksPool.getTask(taskData))
	{
		auto bin = matrixReader->GetBin(taskData.binId);

		uint64_t dataIdx = binsIndicesForCorrection[taskData.binId];
		const uint64_t dataIdxEnd = binsIndicesForCorrection[taskData.binId + 1];

		refresh::statistical_test statistics;
		refresh::scorers scorer;

		std::vector<uint64_t> inMatrixEntry;
		std::ptrdiff_t num_samples = static_cast<std::ptrdiff_t>(params.mkmcParams.samples.size());

		inMatrixEntry.resize(num_samples);

		ProgressBarUpdater progress_bar_updater(*progress_bar, (std::max)(1ull, progress_bar->GetTotal() / 100ull));

		auto kmer_len = params.stage1Params.GetKmerLen();

		kmcdb::DispatchKmerSize<MAX_K>(kmer_len, [&](auto SIZE) {
			kmcdb::CKmer<SIZE> kmer;
			while (bin->NextKmer(kmer, inMatrixEntry.data()))
			{
				size_t outPValuesDataIdx = 0;

				if (statisticsToGeneration.tTest)
				{
					const double tTestPValue = statistics.t_test_n(
						inMatrixEntry.begin(),
						differentialAnalysisPhenotype.begin(),
						num_samples).p_value;

					pValuesData[outPValuesDataIdx++][dataIdx] = tTestPValue;
				}
				if (statisticsToGeneration.wilcoxonRankSum)
				{
					const double wilcoxonRankSumPValue = statistics.mann_whitney_U_test_n(
						inMatrixEntry.begin(),
						differentialAnalysisPhenotype.begin(),
						num_samples).p_value;

					pValuesData[outPValuesDataIdx++][dataIdx] = wilcoxonRankSumPValue;
				}
				if (statisticsToGeneration.anova)
				{
					const double anovaPValue = scorer.anova_n(
						inMatrixEntry.begin(),
						differentialAnalysisPhenotype.begin(),
						differentialAnalysisNClasses,
						num_samples).p_value;

					pValuesData[outPValuesDataIdx][dataIdx] = anovaPValue;
				}

				++dataIdx;
				++progress_bar_updater;
			}
			});

		assert(dataIdx == dataIdxEnd);
	}
}



void StatisticsGenerator::correctPValuesEntries()
{
	typedef StatisticsParams::DifferentialAnalysisCorrectionMethod CorrectionMethod;
	refresh::p_val_correction correction;

	CorrectTaskData correctTaskData;
	while (correctTasksPool.getTask(correctTaskData))
	{

		if (params.statisticsParams.classificationPValueCorrection == CorrectionMethod::Bonferroni)
		{
			correction.bonferroni_n(pValuesData[correctTaskData.algIdx].begin(), pValuesCorrectedData[correctTaskData.algIdx].begin(), pValuesData[correctTaskData.algIdx].size());
		}
		else if (params.statisticsParams.classificationPValueCorrection == CorrectionMethod::HolmBonferroni)
		{
			correction.holm_bonferroni_n(pValuesData[correctTaskData.algIdx].begin(), pValuesCorrectedData[correctTaskData.algIdx].begin(), pValuesData[correctTaskData.algIdx].size());
		}
		else if (params.statisticsParams.classificationPValueCorrection == CorrectionMethod::BenjaminiHochberg)
		{
			correction.benjamini_hochberg_n(pValuesData[correctTaskData.algIdx].begin(), pValuesCorrectedData[correctTaskData.algIdx].begin(), pValuesData[correctTaskData.algIdx].size());
		}
		else if (params.statisticsParams.classificationPValueCorrection == CorrectionMethod::BenjaminiYekutieli)
		{
			correction.benjamini_yekutieli_n(pValuesData[correctTaskData.algIdx].begin(), pValuesCorrectedData[correctTaskData.algIdx].begin(), pValuesData[correctTaskData.algIdx].size());
		}
		else
			assert(false);
	}

}



void StatisticsGenerator::processEntriesAfterCorrection()
{
	TaskData taskData;
	while (tasksPool.getTask(taskData))
	{
		auto bin = matrixReader->GetBin(taskData.binId);

		uint64_t dataIdx = binsIndicesForCorrection[taskData.binId];
		const uint64_t dataIdxEnd = binsIndicesForCorrection[taskData.binId + 1];

		std::unique_ptr<OutputBuffer> normOutputBuffer;

		refresh::normalization_work<uint64_t, double> normalization;
		if (params.statisticsParams.generateNormalization)
		{
			normalization.register_method(params.statisticsParams.normalizationMethod);
			normalization.set_no_series(params.mkmcParams.samples.size());
			normalization.deserialize(params.statisticsParams.normalizationMethod, normalizationData);

			normalization.initialize();

			normOutputBuffer = std::make_unique<OutputBuffer>(*normWriter, getMaxNormLineLength());
		}

		refresh::correlation correlation;
		refresh::statistics_entropy entropyObj;
		refresh::statistical_test statistics;
		refresh::scorers scorer;

		std::vector<uint64_t> inMatrixEntry;
		std::vector<out_kmcdb_value_type> outEntry; // normalized counts and statistics
		std::ptrdiff_t num_samples = static_cast<std::ptrdiff_t>(params.mkmcParams.samples.size());

		inMatrixEntry.resize(num_samples);
		outEntry.resize(statisticsToGeneration.nResults);

		ProgressBarUpdater progress_bar_updater(*progress_bar, (std::max)(1ull, progress_bar->GetTotal() / 100ull));

		auto kmer_len = params.stage1Params.GetKmerLen();
		std::string kmerSequence(kmer_len, ' ');

		kmcdb::DispatchKmerSize<MAX_K>(kmer_len, [&](auto SIZE) {
			kmcdb::CKmer<SIZE> kmer;
			std::unique_ptr<WritingGathererBin<out_kmcdb_value_type>> outBin = gatherer.getBin(taskData.binId);
			while (bin->NextKmer(kmer, inMatrixEntry.data()))
			{
				kmer.to_string(kmer_len, kmerSequence.data());

				if (statisticsToGeneration.normalize)
				{
					normalization.norm_entry(params.statisticsParams.normalizationMethod, inMatrixEntry, outEntry);

					normOutputBuffer->StoreKmer(kmerSequence, outEntry, StoreMethods::AsMatrixRow);

					outEntry.resize(statisticsToGeneration.nResults); // space for statistics
				}

				size_t outStatsEntryIdx = statisticsToGeneration.nResults - statisticsToGeneration.nStatistics;
				size_t outPValuesAlgIdx = 0;

				if (statisticsToGeneration.pearson)
				{
					assert(statisticsToGeneration.normalize);
					const double pearson = refresh::correlation::pearson_n(
						outEntry.begin(),
						correlationPhenotype.begin(),
						num_samples);

					outEntry[outStatsEntryIdx++] = pearson;
				}
				if (statisticsToGeneration.spearman)
				{
					assert(statisticsToGeneration.normalize);
					const double spearman = correlation.spearman_n(
						outEntry.begin(),
						correlationPhenotype.begin(),
						num_samples);

					outEntry[outStatsEntryIdx++] = spearman;
				}
				if (statisticsToGeneration.kendall)
				{
					assert(statisticsToGeneration.normalize);
					const double kendall = refresh::correlation::kendall_tau_n(
						outEntry.begin(),
						correlationPhenotype.begin(),
						num_samples);

					outEntry[outStatsEntryIdx++] = kendall;
				}

				if (statisticsToGeneration.entropy)
				{
					const double entropy = entropyObj.entropy_n(
						inMatrixEntry.begin(),
						num_samples);

					outEntry[outStatsEntryIdx++] = entropy;
				}
				if (statisticsToGeneration.differentialAnalysis)
				{
					if (statisticsToGeneration.tTest)
					{
						outEntry[outStatsEntryIdx++] = pValuesCorrectedData[outPValuesAlgIdx++][dataIdx];
					}
					if (statisticsToGeneration.snr)
					{
						const double snr = statistics.SNR_test_n(
							inMatrixEntry.begin(),
							differentialAnalysisPhenotype.begin(),
							num_samples);

						outEntry[outStatsEntryIdx++] = snr;
					}
					if (statisticsToGeneration.wilcoxonRankSum)
					{
						outEntry[outStatsEntryIdx++] = pValuesCorrectedData[outPValuesAlgIdx++][dataIdx];
					}
					if (statisticsToGeneration.dids)
					{
						const double dids = scorer.dids_n(
							inMatrixEntry.begin(),
							differentialAnalysisPhenotype.begin(),
							differentialAnalysisNClasses,
							num_samples);

						outEntry[outStatsEntryIdx++] = dids;
					}
					if (statisticsToGeneration.anova)
					{
						outEntry[outStatsEntryIdx++] = pValuesCorrectedData[outPValuesAlgIdx++][dataIdx];
					}
				}

				++dataIdx;
				++progress_bar_updater;

				outBin->writeKmer(outEntry, kmer, kmerSequence, inMatrixEntry);
			}
			});
		assert(dataIdx == dataIdxEnd);
	}
}



void StatisticsGenerator::generateStatisticsParallel()
{
	openReaders();
	fillTaskData();

	std::vector<std::string> sample_names;

	std::vector<std::string> cnt_matrix_output_header;
	matrixReader->GetSampleNames(cnt_matrix_output_header);
	assert(!cnt_matrix_output_header.empty());
	if (statisticsToGeneration.normalize)
	{
		sample_names = cnt_matrix_output_header;

		normWriter = std::make_unique<DumpWriter>(params.mkmcParams.outputFileNorm, params.mkmcParams.nThreads > 1);
		normWriter->StoreHeader(sample_names);
	} // otherwise: no normalization in output

	gatherer.initWriting(matrixMetadataReader, sample_names, cnt_matrix_output_header);

	if (params.statisticsParams.generateNormalization)
	{
		MatrixStatsReader stats_reader(params.mkmcParams.normStatsBinFile);
		bool success = false;
		if (params.statisticsParams.normalizationMethod == StatisticsParams::NormalizationMethod::deseq2)
			success = stats_reader.Get(params.statisticsParams.normDeseq2StreamName, normalizationData);
		else if (params.statisticsParams.normalizationMethod == StatisticsParams::NormalizationMethod::frequency_count)
			success = stats_reader.Get(params.statisticsParams.normFrequencyStreamName, normalizationData);
		else if (params.statisticsParams.normalizationMethod == StatisticsParams::NormalizationMethod::quantile)
			success = stats_reader.Get(params.statisticsParams.normQuantileStreamName, normalizationData);

		if (!success)
		{
			std::cerr << "Error: cannot read normalization data\n";
			exit(1);
		}
	}

	if (params.statisticsParams.correctPvalues)
	{
		pValuesData.resize(statisticsToGeneration.nStatisticsWithPValues, std::vector<out_kmcdb_value_type>(binsIndicesForCorrection.back()));

		std::vector<std::thread> threads(params.mkmcParams.nThreads);
		for (uint32_t i_thred = 0; i_thred < params.mkmcParams.nThreads; ++i_thred)
		{
			threads[i_thred] = std::thread([this] { this->gatherPValuesEntriesToCorrection(); });
		}
		for (std::thread& thread : threads)
		{
			thread.join();
		}

		pValuesCorrectedData.resize(statisticsToGeneration.nStatisticsWithPValues, std::vector<out_kmcdb_value_type>(binsIndicesForCorrection.back()));
		for (uint32_t i_thred = 0; i_thred < params.mkmcParams.nThreads; ++i_thred) // probably some threads will be idle
		{
			threads[i_thred] = std::thread([this,i_thred] { this->correctPValuesEntries(); });
		}
		for (std::thread& thread : threads)
		{
			thread.join();
		}

		tasksPool.reset();
		openReaders(); // reopen
		for (uint32_t i_thred = 0; i_thred < params.mkmcParams.nThreads; ++i_thred)
		{
			threads[i_thred] = std::thread([this] { this->processEntriesAfterCorrection(); });
		}
		for (std::thread& thread : threads)
		{
			thread.join();
		}
	}
	else
	{
		std::vector<std::thread> threads(params.mkmcParams.nThreads);
		for (uint32_t i_thred = 0; i_thred < params.mkmcParams.nThreads; ++i_thred)
		{
			threads[i_thred] = std::thread([this] { this->processEntries(); });
		}
		for (std::thread& thread : threads)
		{
			thread.join();
		}
	}
}
