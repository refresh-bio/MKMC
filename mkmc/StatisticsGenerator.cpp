#include "StatisticsGenerator.h"
#include "MatrixStats.h"
#include "DumpWriter.h"
#include <algorithm>



void StatisticsGenerator::fillTaskData()
{
	try
	{
		matrixMetadataReader = std::make_unique<kmcdb::MetadataReader>(params.mkmcParams.outputFilesTemplate + ".kmcdb", false);
		matrixReader = std::make_unique<kmcdb::ReaderSortedPlainForListing<uint64_t>>(*matrixMetadataReader);
	}
	catch (const std::runtime_error& ex)
	{
		std::cerr << "Error: " << ex.what() << std::endl;
		exit(1);
	}
	tasksData.reserve(params.stage1Params.GetNBins());
	std::vector<uint64_t> nOutputKmersPerBin;
	nOutputKmersPerBin.reserve(params.stage1Params.GetNBins());
	
	for (uint32_t i = 0; i < params.stage1Params.GetNBins(); ++i)
	{
		tasksData.push_back(TaskData{ i });
		nOutputKmersPerBin.push_back(matrixReader->GetBin(i)->GetBinMetadata().total_kmers);
	}

	progress_bar = std::make_unique<ProgressBar>(
		params.mkmcParams.verbosity_level == 0 ? 0 : std::accumulate(nOutputKmersPerBin.begin(), nOutputKmersPerBin.end(), 0ull),
		"Computing statistics",
		std::cerr,
		params.mkmcParams.verbosity_level == 0);

	std::sort(tasksData.begin(), tasksData.end(), [&](const TaskData& a, const TaskData& b) { return nOutputKmersPerBin[a.binId] > nOutputKmersPerBin[b.binId]; });
}



StatisticsGenerator::StatisticsGenerator(Params& params) :
	params(params),
	tasksPool(tasksData),
	gatherer(params, statisticsToGeneration),
	correlationPhenotype(params.phenotypes.correlationPhenotype.getPhenotype()),
	differentialAnalysisPhenotype(params.phenotypes.differentialAnalysisPhenotype.getMappedPhenotype()),
	differentialAnalysisNClasses(params.phenotypes.differentialAnalysisPhenotype.getClassesNumber())
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
}



void StatisticsGenerator::operator()()
{
	TaskData taskData;
	while (tasksPool.getTask(taskData))
	{
		auto bin = matrixReader->GetBin(taskData.binId);

		refresh::normalization_work<uint64_t, double> normalization;
		if (params.statisticsParams.generateNormalization)
		{
			normalization.register_method(params.statisticsParams.normalizationMethod);
			normalization.set_no_series(params.mkmcParams.samples.size());
			normalization.deserialize(params.statisticsParams.normalizationMethod, normalizationData);

			normalization.initialize();
		}

		refresh::correlation correlation;
		refresh::statistics_entropy entropyObj;
		refresh::statistical_test statistics;
		refresh::scorers scorer;

		std::vector<uint64_t> inMatrixEntry;
		std::vector<double> outNormMatrixEntry;
		std::vector<out_kmcdb_value_type> outStatsEntry;
		std::ptrdiff_t num_samples = static_cast<std::ptrdiff_t>(params.mkmcParams.samples.size());

		inMatrixEntry.resize(num_samples);
		outNormMatrixEntry.resize(num_samples);
		outStatsEntry.resize(statisticsToGeneration.nStatistics);

		ProgressBarUpdater progress_bar_updater(*progress_bar, (std::max)(1ull, progress_bar->GetTotal() / 100ull));

		auto kmer_len = params.stage1Params.GetKmerLen();
		std::string kmerSequence(kmer_len, ' ');

		kmcdb::DispatchKmerSize<MAX_K>(kmer_len, [&](auto SIZE) {
			kmcdb::CKmer<SIZE> kmer;
			std::unique_ptr<WritingGathererBin<out_kmcdb_value_type>> outBin = gatherer.getBin(taskData.binId);
			while (bin->NextKmer(kmer, inMatrixEntry.data()))
			{
				kmer.to_string(kmer_len, kmerSequence.data());
				normalization.norm_entry(params.statisticsParams.normalizationMethod, inMatrixEntry, outNormMatrixEntry);

				size_t outStatsEntryIdx = 0;

				if (statisticsToGeneration.pearson)
				{
					const double pearson = refresh::correlation::pearson_n(
						outNormMatrixEntry.begin(),
						correlationPhenotype.begin(),
						num_samples);

					outStatsEntry[outStatsEntryIdx++] = pearson;
				}
				if (statisticsToGeneration.spearman)
				{
					const double spearman = correlation.spearman_n(
						outNormMatrixEntry.begin(),
						correlationPhenotype.begin(),
						num_samples);

					outStatsEntry[outStatsEntryIdx++] = spearman;
				}
				if (statisticsToGeneration.kendall)
				{
					const double kendall = refresh::correlation::kendall_tau_n(
						outNormMatrixEntry.begin(),
						correlationPhenotype.begin(),
						num_samples);

					outStatsEntry[outStatsEntryIdx++] = kendall;
				}

				if (statisticsToGeneration.entropy)
				{
					const double entropy = entropyObj.entropy_n(
						outNormMatrixEntry.begin(),
						num_samples);

					outStatsEntry[outStatsEntryIdx++] = entropy;
				}
				if (statisticsToGeneration.differentialAnalysis)
				{
					if (statisticsToGeneration.tTest)
					{
						const double tTestPValue = statistics.t_test_n(
							outNormMatrixEntry.begin(),
							differentialAnalysisPhenotype.begin(),
							num_samples).p_value;

						outStatsEntry[outStatsEntryIdx++] = tTestPValue;
					}
					if (statisticsToGeneration.snr)
					{
						const double snr = statistics.SNR_test_n(
							outNormMatrixEntry.begin(),
							differentialAnalysisPhenotype.begin(),
							num_samples);

						outStatsEntry[outStatsEntryIdx++] = snr;
					}
					if (statisticsToGeneration.wilcoxonRankSum)
					{
						const double wilcoxonRankSumPValue = statistics.mann_whitney_U_test_n(
							outNormMatrixEntry.begin(),
							differentialAnalysisPhenotype.begin(),
							num_samples).p_value;

						outStatsEntry[outStatsEntryIdx++] = wilcoxonRankSumPValue;
					}
					if (statisticsToGeneration.dids)
					{
						const double dids = scorer.dids_n(
							outNormMatrixEntry.begin(),
							differentialAnalysisPhenotype.begin(),
							differentialAnalysisNClasses,
							num_samples);

						outStatsEntry[outStatsEntryIdx++] = dids;
					}
					if (statisticsToGeneration.anova)
					{
						const double anovaPValue = scorer.anova_n(
							outNormMatrixEntry.begin(),
							differentialAnalysisPhenotype.begin(),
							differentialAnalysisNClasses,
							num_samples).p_value;

						outStatsEntry[outStatsEntryIdx++] = anovaPValue;
					}
				}

				++progress_bar_updater;

				outBin->writeKmer(outNormMatrixEntry, kmer, outStatsEntry);
			}
		});
	}
}


void StatisticsGenerator::generateStatisticsParallel()
{
	fillTaskData();

	std::vector<std::string> sample_names;
	matrixReader->GetSampleNames(sample_names);
	assert(!sample_names.empty());

	gatherer.initWriting(matrixMetadataReader, sample_names);

	if (params.statisticsParams.generateNormalization)
	{
		MatrixStatsReader stats_reader(params.mkmcParams.normStatsBinFile);
		bool success = false;
		if (params.statisticsParams.normalizationMethod == StatisticsParams::NormalizationMethod::frequency_count)
			success = stats_reader.Get(params.statisticsParams.normFrequencyStreamName, normalizationData);
		else if (params.statisticsParams.normalizationMethod == StatisticsParams::NormalizationMethod::quantile)
			success = stats_reader.Get(params.statisticsParams.normQuantileStreamName, normalizationData);

		if (!success)
		{
			std::cerr << "Error: cannot read normalization data\n";
			exit(1);
		}
	}

	std::vector<std::thread> threads(params.mkmcParams.nThreads);
	for (uint32_t i_thred = 0; i_thred < params.mkmcParams.nThreads; ++i_thred)
	{
		threads[i_thred] = std::thread([this] { (*this)(); });
	}

	for (std::thread& thread : threads)
	{
		thread.join();
	}
}
