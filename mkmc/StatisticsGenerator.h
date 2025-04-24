#pragma once

#include <fstream>
#include <vector>
#include <string>
#include <memory>
#include "parameters.h"
#include "TasksPool.h"
#include "refresh/statistics/lib/statistics.h"
#define NOMINMAX
#include "progress_bar.hpp"
#include "kmcdb/kmcdb.h"
#include "TextFileWritingUtilities.h"
#include "StatisticsGatherers.h"
#include "DimensionalityReduction.h"


class StatisticsGenerator
{
	Params& params;

	struct TaskData
	{
		uint32_t binId;

		TaskData() :
			binId(static_cast<uint32_t>(-1))
		{}
		TaskData(uint32_t binId) :
			binId(binId)
		{}
	};
	std::vector<TaskData> tasksData;
	TasksPool<TaskData> tasksPool;

	struct CorrectTaskData
	{
		uint32_t algIdx;

		CorrectTaskData() :
			algIdx(static_cast<uint32_t>(-1))
		{}
		CorrectTaskData(uint32_t algIdx) :
			algIdx(algIdx)
		{}
	};
	std::vector<CorrectTaskData> correctTasksData;
	TasksPool<CorrectTaskData> correctTasksPool;

	std::vector<uint64_t> nOutputKmersPerBin;

	using out_kmcdb_value_type = double;
	using cnt_value_type = uint64_t;

	std::unique_ptr<kmcdb::MetadataReader> matrixMetadataReader;
	std::unique_ptr<kmcdb::ReaderSortedPlainForListing<cnt_value_type>> matrixReader;

	WritingGatherer<out_kmcdb_value_type, cnt_value_type> gatherer;
	std::unique_ptr<TextFileWriter> normWriter;

	std::unique_ptr<ProgressBar> progress_bar;

	void openReaders();
	void fillTaskData();


	std::vector<uint8_t> normalizationData;

	const std::vector<out_kmcdb_value_type>& correlationPhenotype;
	const std::vector<uint32_t>& differentialAnalysisPhenotype;
	size_t differentialAnalysisNClasses;

	std::vector<std::vector<out_kmcdb_value_type>> pValuesToCorrect; // first index: algorithm, second: entries
	std::vector<std::vector<out_kmcdb_value_type>> pValuesCorrected; // first index: algorithm, second: entries
	std::vector<std::vector<out_kmcdb_value_type>> additionalValuesOfCorrectedStats; // first index: algorithm, second: entries

	std::vector<uint64_t> binsOffsets; // (of size no. of bins + 1) contains indices of first entries for every bin; the last element contains number of all the entries

	StatisticsToGeneration statisticsToGeneration;

	template<unsigned SIZE>
	void initKeepNLargest(KeepNLargestCollection<SIZE>& keepNLargestCollection);

	template<unsigned SIZE>
	void processEntries(KeepNLargestCollectionGlobal<SIZE, out_kmcdb_value_type, cnt_value_type>& keepNLargestCollectionGlobal,
		DimensionalityReduction& dimensionalityReduction);

	template<unsigned SIZE>
	void processEntriesWhenCorrection(KeepNLargestCollectionGlobal<SIZE, out_kmcdb_value_type, cnt_value_type>& keepNLargestCollectionGlobal,
		DimensionalityReduction& dimensionalityReduction);

	void correctPValuesEntries();

	template<unsigned SIZE>
	void safeCorrectedPValuesEntries();

	bool readNormalizationData();

public:
	StatisticsGenerator(Params& params);

	void generateStatisticsParallel();
};


template<unsigned SIZE>
void StatisticsGenerator::initKeepNLargest(KeepNLargestCollection<SIZE>& keepNLargestCollection)
{
	if (params.statisticsParams.nTop == 0)
		return;

	using KeepTopNLargestABS_T = typename KeepNLargestCollection<SIZE>::KeepTopNLargestABS_T;
	using KeepTopNLargestPlain_T = typename KeepNLargestCollection<SIZE>::KeepTopNLargestPlain_T;

	if (statisticsToGeneration.pearson)
		keepNLargestCollection.pearson = std::make_unique<KeepTopNLargestABS_T>(params.statisticsParams.nTop);

	if (statisticsToGeneration.spearman)
		keepNLargestCollection.spearman = std::make_unique<KeepTopNLargestABS_T>(params.statisticsParams.nTop);

	if (statisticsToGeneration.kendall)
		keepNLargestCollection.kendall = std::make_unique<KeepTopNLargestABS_T>(params.statisticsParams.nTop);

	if (statisticsToGeneration.entropy)
		keepNLargestCollection.entropy = std::make_unique<KeepTopNLargestPlain_T>(params.statisticsParams.nTop);

	if (statisticsToGeneration.snr)
		keepNLargestCollection.snr = std::make_unique<KeepTopNLargestPlain_T>(params.statisticsParams.nTop);

	if (statisticsToGeneration.unnormalizedSnr)
		keepNLargestCollection.unnormalizedSnr = std::make_unique<KeepTopNLargestPlain_T>(params.statisticsParams.nTop);

	if (statisticsToGeneration.dids)
		keepNLargestCollection.dids = std::make_unique<KeepTopNLargestPlain_T>(params.statisticsParams.nTop);

}

template<unsigned SIZE>
void StatisticsGenerator::processEntries(KeepNLargestCollectionGlobal<SIZE, out_kmcdb_value_type, cnt_value_type>& keepNLargestCollectionGlobal,
	DimensionalityReduction& dimensionalityReduction)
{
	KeepNLargestCollection<SIZE> keepNLargestCollection;

	initKeepNLargest(keepNLargestCollection);

	const std::size_t num_samples = params.mkmcParams.samples.size();

	std::vector<cnt_value_type> inMatrixEntry;
	std::vector<out_kmcdb_value_type> outNormEntry; // normalized stats
	std::vector<out_kmcdb_value_type> outStatsEntry; // statistics
	inMatrixEntry.resize(num_samples);
	outNormEntry.resize(num_samples);
	outStatsEntry.resize(statisticsToGeneration.nStatistics + statisticsToGeneration.nAdditionalValuesOfCorrectedStats);

	const auto first_col_len = params.stage1Params.GetKmerLen();
	std::string kmerSequence(first_col_len, ' ');

	TaskData taskData;
	while (tasksPool.getTask(taskData))
	{
		auto bin = matrixReader->GetBin(taskData.binId);

		std::unique_ptr<MatrixOutputBuffer<out_kmcdb_value_type>> normOutputBuffer;

		refresh::normalization_work<cnt_value_type, out_kmcdb_value_type> normalization;
		if (params.statisticsParams.generateNormalization)
		{
			normalization.register_method(params.statisticsParams.normalizationMethod);
			normalization.set_no_series(params.mkmcParams.samples.size());
			normalization.deserialize(params.statisticsParams.normalizationMethod, normalizationData);

			normalization.initialize();

			normOutputBuffer = std::make_unique<MatrixOutputBuffer<out_kmcdb_value_type>>(*normWriter, params.stage1Params.GetKmerLen(), num_samples);
		}

		refresh::correlation correlation;
		refresh::statistics_entropy entropyObj;
		refresh::statistical_test statistics;
		refresh::scorers scorer;

		ProgressBarUpdater progress_bar_updater(*progress_bar, (std::max)(1ull, progress_bar->GetTotal() / 100ull));

		kmcdb::CKmer<SIZE> kmer;
		std::unique_ptr<WritingGathererBin<out_kmcdb_value_type, cnt_value_type>> outGahtererBin = gatherer.getBin(taskData.binId);

		std::vector<double> inMatrixEntryScaled(num_samples); // for t-test
		for (auto kmer_idx = binsOffsets[taskData.binId]; bin->NextKmer(kmer, inMatrixEntry.data()); ++kmer_idx)
		{
			kmer.to_string(first_col_len, kmerSequence.data());

			if (statisticsToGeneration.normalize)
			{
				normalization.norm_entry(params.statisticsParams.normalizationMethod, inMatrixEntry, outNormEntry);

				dimensionalityReduction.add(kmer_idx, outNormEntry);
				normOutputBuffer->StoreKmer(kmerSequence, outNormEntry);
			}

			size_t outStatsEntryIdx = 0;

			if (statisticsToGeneration.pearson)
			{
				assert(statisticsToGeneration.normalize);
				const double pearson = refresh::correlation::pearson_n(
					outNormEntry.begin(),
					correlationPhenotype.begin(),
					num_samples);

				outStatsEntry[outStatsEntryIdx++] = pearson;
			}
			if (statisticsToGeneration.spearman)
			{
				assert(statisticsToGeneration.normalize);
				const double spearman = correlation.spearman_n(
					outNormEntry.begin(),
					correlationPhenotype.begin(),
					num_samples);

				outStatsEntry[outStatsEntryIdx++] = spearman;
			}
			if (statisticsToGeneration.kendall)
			{
				assert(statisticsToGeneration.normalize);
				const double kendall = refresh::correlation::kendall_tau_n(
					outNormEntry.begin(),
					correlationPhenotype.begin(),
					num_samples);

				outStatsEntry[outStatsEntryIdx++] = kendall;
			}

			if (statisticsToGeneration.entropy)
			{
				const double entropy = entropyObj.entropy_scaled_n(
					inMatrixEntry.begin(),
					num_samples);

				outStatsEntry[outStatsEntryIdx++] = entropy;
			}
			if (statisticsToGeneration.differentialAnalysis)
			{
				if (statisticsToGeneration.tTest)
				{
					for (size_t it = 0; it < num_samples; ++it)
						inMatrixEntryScaled[it] = std::log2(inMatrixEntry[it] + 1);

					const auto tTest = statistics.t_test_n(
						inMatrixEntryScaled.begin(),
						differentialAnalysisPhenotype.begin(),
						num_samples,
						false);

					outStatsEntry[outStatsEntryIdx++] = tTest.p_value;
					outStatsEntry[outStatsEntryIdx++] = tTest.df;
					outStatsEntry[outStatsEntryIdx++] = tTest.statistic;
				}
				if (statisticsToGeneration.snr)
				{
					const double snr = statistics.SNR_test_n(
						outNormEntry.begin(),
						differentialAnalysisPhenotype.begin(),
						num_samples);

					outStatsEntry[outStatsEntryIdx++] = snr;
				}
				if (statisticsToGeneration.unnormalizedSnr)
				{
					const double snr = statistics.SNR_test_n(
						inMatrixEntry.begin(),
						differentialAnalysisPhenotype.begin(),
						num_samples);

					outStatsEntry[outStatsEntryIdx++] = snr;
				}
				if (statisticsToGeneration.wrs)
				{
					const auto wrs = statistics.mann_whitney_U_test_n(
						outNormEntry.begin(),
						differentialAnalysisPhenotype.begin(),
						num_samples);

					outStatsEntry[outStatsEntryIdx++] = wrs.p_value;
					outStatsEntry[outStatsEntryIdx++] = wrs.statistic_U1;
					outStatsEntry[outStatsEntryIdx++] = wrs.statistic_U2;
				}
				if (statisticsToGeneration.dids)
				{
					const double dids = scorer.dids_n(
						outNormEntry.begin(),
						differentialAnalysisPhenotype.begin(),
						differentialAnalysisNClasses,
						num_samples);

					outStatsEntry[outStatsEntryIdx++] = dids;
				}
				if (statisticsToGeneration.anova)
				{
					const auto anova = scorer.anova_n(
						outNormEntry.begin(),
						differentialAnalysisPhenotype.begin(),
						differentialAnalysisNClasses,
						num_samples);

					outStatsEntry[outStatsEntryIdx++] = anova.p_value;
					outStatsEntry[outStatsEntryIdx++] = anova.statistic;
				}
			}
			assert(outStatsEntryIdx == statisticsToGeneration.nStatistics + statisticsToGeneration.nAdditionalValuesOfCorrectedStats);

			++progress_bar_updater;

			outGahtererBin->writeKmer(outStatsEntry, kmer, kmerSequence, inMatrixEntry, &keepNLargestCollection);
		}
	}
	keepNLargestCollectionGlobal.Add(keepNLargestCollection);
}



template<unsigned SIZE>
void StatisticsGenerator::processEntriesWhenCorrection(KeepNLargestCollectionGlobal<SIZE, out_kmcdb_value_type, cnt_value_type>& keepNLargestCollectionGlobal,
	DimensionalityReduction& dimensionalityReduction)
{
	assert(statisticsToGeneration.differentialAnalysis);
	assert(params.statisticsParams.generateNormalization);

	KeepNLargestCollection<SIZE> keepNLargestCollection;

	initKeepNLargest(keepNLargestCollection);

	const std::size_t num_samples = params.mkmcParams.samples.size();

	std::vector<cnt_value_type> inMatrixEntry;
	std::vector<out_kmcdb_value_type> outNormEntry; // normalized counts
	std::vector<out_kmcdb_value_type> outStatsEntry; // statistics
	inMatrixEntry.resize(num_samples);
	outNormEntry.resize(num_samples);
	outStatsEntry.resize(statisticsToGeneration.nStatistics - statisticsToGeneration.nStatisticsWithPValues);

	const auto first_col_len = params.stage1Params.GetKmerLen();
	std::string kmerSequence(first_col_len, ' ');

	TaskData taskData;
	while (tasksPool.getTask(taskData))
	{
		auto bin = matrixReader->GetBin(taskData.binId);

		uint64_t outPValuesToCorrectIdx = binsOffsets[taskData.binId];
		const uint64_t outPValuesToCorrectIdxEnd = binsOffsets[taskData.binId + 1];

		std::unique_ptr<MatrixOutputBuffer<out_kmcdb_value_type>> normOutputBuffer = std::make_unique<MatrixOutputBuffer<out_kmcdb_value_type>>(*normWriter, params.stage1Params.GetKmerLen(), num_samples);

		assert(params.statisticsParams.generateNormalization);
		refresh::normalization_work<cnt_value_type, out_kmcdb_value_type> normalization;
		normalization.register_method(params.statisticsParams.normalizationMethod);
		normalization.set_no_series(params.mkmcParams.samples.size());
		normalization.deserialize(params.statisticsParams.normalizationMethod, normalizationData);

		normalization.initialize();

		refresh::correlation correlation;
		refresh::statistics_entropy entropyObj;
		refresh::statistical_test statistics;
		refresh::scorers scorer;

		ProgressBarUpdater progress_bar_updater(*progress_bar, (std::max)(1ull, progress_bar->GetTotal() / 100ull));

		kmcdb::CKmer<SIZE> kmer;
		std::unique_ptr<WritingGathererBin<out_kmcdb_value_type, cnt_value_type>> outGathererBin = gatherer.getBin(taskData.binId, false, true);

		std::vector<double> inMatrixEntryScaled(num_samples); // for t-test
		while (bin->NextKmer(kmer, inMatrixEntry.data()))
		{
			kmer.to_string(first_col_len, kmerSequence.data());

			normalization.norm_entry(params.statisticsParams.normalizationMethod, inMatrixEntry, outNormEntry);

			dimensionalityReduction.add(outPValuesToCorrectIdx, outNormEntry);
			normOutputBuffer->StoreKmer(kmerSequence, outNormEntry);

			size_t outStatsIdx = 0;
			size_t outPValuesToCorrectAlg = 0;
			size_t outAdditionalValuesIdx = 0;

			if (statisticsToGeneration.pearson)
			{
				assert(statisticsToGeneration.normalize);
				const double pearson = refresh::correlation::pearson_n(
					outNormEntry.begin(),
					correlationPhenotype.begin(),
					num_samples);

				outStatsEntry[outStatsIdx++] = pearson;
			}
			if (statisticsToGeneration.spearman)
			{
				assert(statisticsToGeneration.normalize);
				const double spearman = correlation.spearman_n(
					outNormEntry.begin(),
					correlationPhenotype.begin(),
					num_samples);

				outStatsEntry[outStatsIdx++] = spearman;
			}
			if (statisticsToGeneration.kendall)
			{
				assert(statisticsToGeneration.normalize);
				const double kendall = refresh::correlation::kendall_tau_n(
					outNormEntry.begin(),
					correlationPhenotype.begin(),
					num_samples);

				outStatsEntry[outStatsIdx++] = kendall;
			}

			if (statisticsToGeneration.entropy)
			{
				const double entropy = entropyObj.entropy_scaled_n(
					inMatrixEntry.begin(),
					num_samples);

				outStatsEntry[outStatsIdx++] = entropy;
			}

			if (statisticsToGeneration.tTest)
			{
				for (size_t it = 0; it < num_samples; ++it)
					inMatrixEntryScaled[it] = std::log2(inMatrixEntry[it] + 1);

				const auto tTest = statistics.t_test_n(
					inMatrixEntryScaled.begin(),
					differentialAnalysisPhenotype.begin(),
					num_samples,
					false);

				pValuesToCorrect[outPValuesToCorrectAlg++][outPValuesToCorrectIdx] = tTest.p_value;
				additionalValuesOfCorrectedStats[outAdditionalValuesIdx++][outPValuesToCorrectIdx] = tTest.df;
				additionalValuesOfCorrectedStats[outAdditionalValuesIdx++][outPValuesToCorrectIdx] = tTest.statistic;
			}
			if (statisticsToGeneration.snr)
			{
				const double snr = statistics.SNR_test_n(
					outNormEntry.begin(),
					differentialAnalysisPhenotype.begin(),
					num_samples);

				outStatsEntry[outStatsIdx++] = snr;
			}
			if (statisticsToGeneration.unnormalizedSnr)
			{
				const double snr = statistics.SNR_test_n(
					inMatrixEntry.begin(),
					differentialAnalysisPhenotype.begin(),
					num_samples);

				outStatsEntry[outStatsIdx++] = snr;
			}
			if (statisticsToGeneration.wrs)
			{
				const auto wrs = statistics.mann_whitney_U_test_n(
					outNormEntry.begin(),
					differentialAnalysisPhenotype.begin(),
					num_samples);

				pValuesToCorrect[outPValuesToCorrectAlg++][outPValuesToCorrectIdx] = wrs.p_value;
				additionalValuesOfCorrectedStats[outAdditionalValuesIdx++][outPValuesToCorrectIdx] = wrs.statistic_U1;
				additionalValuesOfCorrectedStats[outAdditionalValuesIdx++][outPValuesToCorrectIdx] = wrs.statistic_U2;
			}
			if (statisticsToGeneration.dids)
			{
				const double dids = scorer.dids_n(
					outNormEntry.begin(),
					differentialAnalysisPhenotype.begin(),
					differentialAnalysisNClasses,
					num_samples);

				outStatsEntry[outStatsIdx++] = dids;
			}
			if (statisticsToGeneration.anova)
			{
				const auto anova = scorer.anova_n(
					outNormEntry.begin(),
					differentialAnalysisPhenotype.begin(),
					differentialAnalysisNClasses,
					num_samples);

				pValuesToCorrect[outPValuesToCorrectAlg++][outPValuesToCorrectIdx] = anova.p_value;
				additionalValuesOfCorrectedStats[outAdditionalValuesIdx++][outPValuesToCorrectIdx] = anova.statistic;
			}
			assert(outPValuesToCorrectAlg == statisticsToGeneration.nStatisticsWithPValues);

			++outPValuesToCorrectIdx;
			++progress_bar_updater;

			outGathererBin->writeKmer(outStatsEntry, kmer, kmerSequence, inMatrixEntry, &keepNLargestCollection);
		}

		assert(outPValuesToCorrectIdx == outPValuesToCorrectIdxEnd);
	}
	keepNLargestCollectionGlobal.Add(keepNLargestCollection);
}



template<unsigned SIZE>
void StatisticsGenerator::safeCorrectedPValuesEntries()
{
	const std::size_t num_samples = params.mkmcParams.samples.size();

	std::vector<cnt_value_type> inMatrixEntry;
	std::vector<out_kmcdb_value_type> outStatsEntry; // statistics
	inMatrixEntry.resize(num_samples);
	outStatsEntry.resize(statisticsToGeneration.nStatistics + statisticsToGeneration.nAdditionalValuesOfCorrectedStats);

	const auto first_col_len = params.stage1Params.GetKmerLen();
	std::string kmerSequence(first_col_len, ' ');

	TaskData taskData;
	while (tasksPool.getTask(taskData))
	{
		auto bin = matrixReader->GetBin(taskData.binId);

		uint64_t outPValuesToCorrectIdx = binsOffsets[taskData.binId];
		const uint64_t outPValuesToCorrectIdxEnd = binsOffsets[taskData.binId + 1];

		ProgressBarUpdater progress_bar_updater(*progress_bar, (std::max)(1ull, progress_bar->GetTotal() / 100ull));

		kmcdb::CKmer<SIZE> kmer;
		std::unique_ptr<WritingGathererBin<out_kmcdb_value_type, cnt_value_type>> outGathererBin = gatherer.getBin(taskData.binId, true, false);
		while (bin->NextKmer(kmer, inMatrixEntry.data()))
		{
			kmer.to_string(first_col_len, kmerSequence.data());

			size_t outPvaluesIdx = 0;
			size_t outAdditionalValuesIdx = 0;
			if (statisticsToGeneration.differentialAnalysis)
			{
				if (statisticsToGeneration.tTest)
				{
					outStatsEntry[outPvaluesIdx + outAdditionalValuesIdx] = pValuesCorrected[outPvaluesIdx][outPValuesToCorrectIdx];
					++outPvaluesIdx;

					// additional values
					outStatsEntry[outPvaluesIdx + outAdditionalValuesIdx] = additionalValuesOfCorrectedStats[outAdditionalValuesIdx][outPValuesToCorrectIdx];
					++outAdditionalValuesIdx;
					outStatsEntry[outPvaluesIdx + outAdditionalValuesIdx] = additionalValuesOfCorrectedStats[outAdditionalValuesIdx][outPValuesToCorrectIdx];
					++outAdditionalValuesIdx;
				}
				if (statisticsToGeneration.wrs)
				{
					outStatsEntry[outPvaluesIdx + outAdditionalValuesIdx] = pValuesCorrected[outPvaluesIdx][outPValuesToCorrectIdx];
					++outPvaluesIdx;

					// additional values
					outStatsEntry[outPvaluesIdx + outAdditionalValuesIdx] = additionalValuesOfCorrectedStats[outAdditionalValuesIdx][outPValuesToCorrectIdx];
					++outAdditionalValuesIdx;
					outStatsEntry[outPvaluesIdx + outAdditionalValuesIdx] = additionalValuesOfCorrectedStats[outAdditionalValuesIdx][outPValuesToCorrectIdx];
					++outAdditionalValuesIdx;
				}
				if (statisticsToGeneration.anova)
				{
					outStatsEntry[outPvaluesIdx + outAdditionalValuesIdx] = pValuesCorrected[outPvaluesIdx][outPValuesToCorrectIdx];
					++outPvaluesIdx;

					// additional value
					outStatsEntry[outPvaluesIdx + outAdditionalValuesIdx] = additionalValuesOfCorrectedStats[outAdditionalValuesIdx][outPValuesToCorrectIdx];
					++outAdditionalValuesIdx;
				}
			}
			assert(outPvaluesIdx + outAdditionalValuesIdx == statisticsToGeneration.nStatisticsWithPValues + statisticsToGeneration.nAdditionalValuesOfCorrectedStats);

			++outPValuesToCorrectIdx;
			++progress_bar_updater;

			outGathererBin->writeKmer(outStatsEntry, kmer, kmerSequence, inMatrixEntry);
		}
		assert(outPValuesToCorrectIdx == outPValuesToCorrectIdxEnd);
	}
}
