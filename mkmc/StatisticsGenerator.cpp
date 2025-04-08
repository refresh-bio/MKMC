#include "StatisticsGenerator.h"
#include "MatrixStats.h"
#include "DumpWriter.h"
#include "DimensionalityReduction.h"
#include "Deseq2Learner.h"
#include <algorithm>



void StatisticsGenerator::openReaders()
{
	try
	{
		matrixMetadataReader = std::make_unique<kmcdb::MetadataReader>(params.mkmcParams.outputMatrixBinFile, false);
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
	tasksData.reserve(matrixMetadataReader->GetConfig().num_bins);
	if (params.statisticsParams.correctPvalues)
		correctTasksData.reserve(statisticsToGeneration.nStatisticsWithPValues);

	nOutputKmersPerBin.reserve(matrixMetadataReader->GetConfig().num_bins);
	
	for (uint32_t i = 0; i < matrixMetadataReader->GetConfig().num_bins; ++i)
	{
		tasksData.push_back(TaskData{ i });
		nOutputKmersPerBin.push_back(matrixReader->GetBin(i)->GetBinMetadata().total_kmers);
	}

	if (params.statisticsParams.correctPvalues)
		for (uint32_t i = 0; i < statisticsToGeneration.nStatisticsWithPValues; ++i)
			correctTasksData.push_back(CorrectTaskData{ i });

	binsOffsets.resize(matrixMetadataReader->GetConfig().num_bins + 1);

	for (size_t it = 1; it < binsOffsets.size(); ++it) // fragmentsBegins[0] = 0
	{
		binsOffsets[it] = binsOffsets[it - 1] + nOutputKmersPerBin[it - 1]; // increase previous index by a size of the next bin
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
	// correctTasksData does not need to be sorted
}



bool StatisticsGenerator::readNormalizationData()
{
	MatrixStatsReader stats_reader(params.mkmcParams.normLearningBinFile);
	bool success = false;
	if (params.statisticsParams.normalizationMethod == StatisticsParams::NormalizationMethod::deseq2) {
		success = stats_reader.Get(params.statisticsParams.normDeseq2StreamName, normalizationData);
		if (!success) {
			// try to open file supplemented with DESeq2
			try { // will be useful after modularization
				MatrixStatsReader stats_reader_supplemented(params.mkmcParams.normLearningBinFileSupplemented);
				success = stats_reader_supplemented.Get(params.statisticsParams.normDeseq2StreamName, normalizationData);
			}
			catch (...) {
				// do nothing, because missing file is not a problem symptom
			}
			if (success) {
				std::cerr << "Info: previously supplemented learning data for DESeq2 properly opened\n";
			}
			else { // learn also for DESeq2, if not learned eariler; will be useful after modularization
				std::cerr << "Info: DESeq2 learning data is missing; it will be supplemented\n";
				params.statisticsParams.normalizationLearningWasSupplemented = true;

				Deseq2LearnerRunner deseq2LearnerRunner(params);
				DispatchKmerSize(params.stage1Params.GetKmerLen(), deseq2LearnerRunner);

				// try to open file lately supplemented with DESeq2
				try {
					MatrixStatsReader stats_reader_currently_supplemented(params.mkmcParams.normLearningBinFileSupplemented);
					success = stats_reader_currently_supplemented.Get(params.statisticsParams.normDeseq2StreamName, normalizationData);
				}
				catch (...) {
					// do nothing, because success == false cause following error message
				}
			}
		}
	}
	else if (params.statisticsParams.normalizationMethod == StatisticsParams::NormalizationMethod::frequency_count)
		success = stats_reader.Get(params.statisticsParams.normFrequencyStreamName, normalizationData);
	else if (params.statisticsParams.normalizationMethod == StatisticsParams::NormalizationMethod::quantile)
		success = stats_reader.Get(params.statisticsParams.normQuantileStreamName, normalizationData);

	return success;
}



StatisticsGenerator::StatisticsGenerator(Params& params) :
	params(params),
	tasksPool(tasksData),
	correctTasksPool(correctTasksData),
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

	statisticsToGeneration.normalize = params.statisticsParams.generateNormalization;

	statisticsToGeneration.pearson = is_correlation_method(StatisticsParams::CorrelationMethod::Pearson);
	statisticsToGeneration.spearman = is_correlation_method(StatisticsParams::CorrelationMethod::Spearman);
	statisticsToGeneration.kendall = is_correlation_method(StatisticsParams::CorrelationMethod::Kendall);

	statisticsToGeneration.entropy = params.statisticsParams.generateEntropy;
	if (statisticsToGeneration.entropy)
		++statisticsToGeneration.nStatistics;

	statisticsToGeneration.tTest = is_differential_analysis_method(StatisticsParams::DifferentialAnalysisMethod::TTest);
	statisticsToGeneration.snr = is_differential_analysis_method(StatisticsParams::DifferentialAnalysisMethod::SNR);
	if (params.mkmcParams.generateForNonNormalized) // as is_differential_analysis_method has side effects, it should be called conditionally here
		statisticsToGeneration.unnormalizedSnr = is_differential_analysis_method(StatisticsParams::DifferentialAnalysisMethod::SNR);
	statisticsToGeneration.wilcoxonRankSum = is_differential_analysis_method(StatisticsParams::DifferentialAnalysisMethod::WilcoxonRankSum);

	statisticsToGeneration.dids = is_differential_analysis_method(StatisticsParams::DifferentialAnalysisMethod::DIDS);
	statisticsToGeneration.anova = is_differential_analysis_method(StatisticsParams::DifferentialAnalysisMethod::ANOVA);

	statisticsToGeneration.differentialAnalysis = statisticsToGeneration.tTest || statisticsToGeneration.snr || statisticsToGeneration.wilcoxonRankSum || statisticsToGeneration.dids || statisticsToGeneration.anova;

	if (statisticsToGeneration.tTest)
	{
		++statisticsToGeneration.nStatisticsWithPValues;
		statisticsToGeneration.nAdditionalValuesOfCorrectedStats += 2;
	}
	if (statisticsToGeneration.wilcoxonRankSum)
	{
		++statisticsToGeneration.nStatisticsWithPValues;
		statisticsToGeneration.nAdditionalValuesOfCorrectedStats += 2;
	}
	if (statisticsToGeneration.anova)
	{
		++statisticsToGeneration.nStatisticsWithPValues;
		++statisticsToGeneration.nAdditionalValuesOfCorrectedStats;
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
			correction.bonferroni_n(pValuesToCorrect[correctTaskData.algIdx].begin(), pValuesCorrected[correctTaskData.algIdx].begin(), pValuesToCorrect[correctTaskData.algIdx].size());
		}
		else if (params.statisticsParams.classificationPValueCorrection == CorrectionMethod::HolmBonferroni)
		{
			correction.holm_bonferroni_n(pValuesToCorrect[correctTaskData.algIdx].begin(), pValuesCorrected[correctTaskData.algIdx].begin(), pValuesToCorrect[correctTaskData.algIdx].size());
		}
		else if (params.statisticsParams.classificationPValueCorrection == CorrectionMethod::BenjaminiHochberg)
		{
			correction.benjamini_hochberg_n(pValuesToCorrect[correctTaskData.algIdx].begin(), pValuesCorrected[correctTaskData.algIdx].begin(), pValuesToCorrect[correctTaskData.algIdx].size());
		}
		else if (params.statisticsParams.classificationPValueCorrection == CorrectionMethod::BenjaminiYekutieli)
		{
			correction.benjamini_yekutieli_n(pValuesToCorrect[correctTaskData.algIdx].begin(), pValuesCorrected[correctTaskData.algIdx].begin(), pValuesToCorrect[correctTaskData.algIdx].size());
		}
		else
			assert(false);
	}

}




void StatisticsGenerator::generateStatisticsParallel()
{
	openReaders();
	fillTaskData();

	std::vector<std::string> samples_names;
	matrixReader->GetSampleNames(samples_names);
	assert(!samples_names.empty());
	if (statisticsToGeneration.normalize)
	{
		normWriter = std::make_unique<DumpWriter>(params.mkmcParams.outputFileNorm, params.mkmcParams.nThreads > 1);
		normWriter->StoreHeader(samples_names);
	} // otherwise: no normalization in output

	gatherer.initWriting(matrixMetadataReader, samples_names);

	if (params.statisticsParams.generateNormalization)
	{
		if (!readNormalizationData())
		{
			std::cerr << "Error: cannot read normalization data\n";
			exit(1);
		}
	}

	if (params.statisticsParams.correctPvalues)
	{
		pValuesToCorrect.resize(statisticsToGeneration.nStatisticsWithPValues, std::vector<out_kmcdb_value_type>(binsOffsets.back()));
		additionalValuesOfCorrectedStats.resize(statisticsToGeneration.nAdditionalValuesOfCorrectedStats, std::vector<out_kmcdb_value_type>(binsOffsets.back()));

		std::vector<std::thread> threads(params.mkmcParams.nThreads);

		kmcdb::DispatchKmerSize<MAX_K>(params.stage1Params.GetKmerLen(), [&](auto SIZE) {
			DimensionalityReduction dimensionalityReduction(params,
				samples_names,
				binsOffsets.back() //number of k-mers
			);
			KeepNLargestCollectionGlobal<SIZE> keepNLargestCollectionGlobal;

			for (uint32_t i_thred = 0; i_thred < params.mkmcParams.nThreads; ++i_thred)
			{
				threads[i_thred] = std::thread([this, &keepNLargestCollectionGlobal, &dimensionalityReduction]
					{ this->processEntriesWhenCorrection<decltype(SIZE)::value>(keepNLargestCollectionGlobal, dimensionalityReduction); });
			}
			for (std::thread& thread : threads)
			{
				thread.join();
			}

			dimensionalityReduction.runAndStore();
			keepNLargestCollectionGlobal.Flush(params, samples_names); // samples names as matrix header
		});

		pValuesCorrected.resize(statisticsToGeneration.nStatisticsWithPValues, std::vector<out_kmcdb_value_type>(binsOffsets.back()));
		for (uint32_t i_thred = 0; i_thred < params.mkmcParams.nThreads; ++i_thred) // probably some threads will be idle
		{
			threads[i_thred] = std::thread([this] { this->correctPValuesEntries(); });
		}
		for (std::thread& thread : threads)
		{
			thread.join();
		}

		tasksPool.reset();
		openReaders(); // reopen
		kmcdb::DispatchKmerSize<MAX_K>(params.stage1Params.GetKmerLen(), [&](auto SIZE) {
			for (uint32_t i_thred = 0; i_thred < params.mkmcParams.nThreads; ++i_thred)
			{
				threads[i_thred] = std::thread([this]
					{ this->safeCorrectedPValuesEntries<decltype(SIZE)::value>(); });
			}
			for (std::thread& thread : threads)
			{
				thread.join();
			}
		});
	}
	else
	{
		kmcdb::DispatchKmerSize<MAX_K>(params.stage1Params.GetKmerLen(), [&](auto SIZE) {
			DimensionalityReduction dimensionalityReduction(params,
				samples_names,
				binsOffsets.back() //number of k-mers
				);
			KeepNLargestCollectionGlobal<SIZE> keepNLargestCollectionGlobal;

			std::vector<std::thread> threads(params.mkmcParams.nThreads);
			for (uint32_t i_thred = 0; i_thred < params.mkmcParams.nThreads; ++i_thred)
			{
				threads[i_thred] = std::thread([this, &keepNLargestCollectionGlobal,&dimensionalityReduction]
					{ this->processEntries<decltype(SIZE)::value>(keepNLargestCollectionGlobal, dimensionalityReduction); });
			}
			for (std::thread& thread : threads)
			{
				thread.join();
			}

			dimensionalityReduction.runAndStore();
			keepNLargestCollectionGlobal.Flush(params, samples_names); // samples names as matrix header
		});
	}
}
