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
#include "DumpWriter.h"
#include "StatisticsGatherers.h"


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

	std::unique_ptr<kmcdb::MetadataReader> matrixMetadataReader;
	std::unique_ptr<kmcdb::ReaderSortedPlainForListing<uint64_t>> matrixReader;

	using out_kmcdb_value_type = double;

	WritingGatherer<out_kmcdb_value_type> gatherer;
	std::unique_ptr<DumpWriter> normWriter;

	std::unique_ptr<ProgressBar> progress_bar;

	void openReaders();
	void fillTaskData();


	std::vector<uint8_t> normalizationData;

	const std::vector<double>& correlationPhenotype;
	const std::vector<uint32_t>& differentialAnalysisPhenotype;
	size_t differentialAnalysisNClasses;

	std::vector<std::vector<out_kmcdb_value_type>> pValuesData; // first index: algorithm, second: entries
	std::vector<std::vector<out_kmcdb_value_type>> pValuesCorrectedData; // first index: algorithm, second: entries

	//mkokot_TODO: I am using this also for umap, so maybe change a name of this variable, like binsOffsets
	std::vector<uint64_t> binsIndicesForCorrection; // (of size no. of bins + 1) contains indices of first entries for every bin; the last element contains number of all the entries

	StatisticsToGeneration statisticsToGeneration;

	size_t getMaxNormLineLength() const
	{
		return params.stage1Params.GetKmerLen() + 1 + params.mkmcParams.samples.size() * (refresh::numeric_conversion_max_length<out_kmcdb_value_type>() + 1);
	}

	template<unsigned SIZE>
	void initKeepNLargest(KeepNLargestCollection<SIZE>& keepNLargestCollection);

	template<unsigned SIZE>
	void processEntries(KeepNLargestCollectionGlobal<SIZE>& keepNLargestCollectionGlobal,
		refresh::umap_direct<double>* umap);

	template<unsigned SIZE>
	void processEntriesWhenCorrection(KeepNLargestCollectionGlobal<SIZE>& keepNLargestCollectionGlobal, refresh::umap_direct<double>* umap);

	void correctPValuesEntries();

	template<unsigned SIZE>
	void safeCorrectedPValuesEntries();

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
		keepNLargestCollection.snr = std::make_unique<KeepTopNLargestABS_T>(params.statisticsParams.nTop);

	if (statisticsToGeneration.dids)
		keepNLargestCollection.dids = std::make_unique<KeepTopNLargestPlain_T>(params.statisticsParams.nTop);

}

template<unsigned SIZE>
void StatisticsGenerator::processEntries(KeepNLargestCollectionGlobal<SIZE>& keepNLargestCollectionGlobal,
	refresh::umap_direct<double>* umap)
{
	KeepNLargestCollection<SIZE> keepNLargestCollection;

	initKeepNLargest(keepNLargestCollection);

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
		std::vector<out_kmcdb_value_type> outNormEntry; // normalized stats
		std::vector<out_kmcdb_value_type> outStatsEntry; // statistics
		std::ptrdiff_t num_samples = static_cast<std::ptrdiff_t>(params.mkmcParams.samples.size());

		inMatrixEntry.resize(num_samples);
		outNormEntry.resize(num_samples);
		outStatsEntry.resize(statisticsToGeneration.nStatistics);

		ProgressBarUpdater progress_bar_updater(*progress_bar, (std::max)(1ull, progress_bar->GetTotal() / 100ull));

		auto kmer_len = params.stage1Params.GetKmerLen();
		std::string kmerSequence(kmer_len, ' ');

		kmcdb::CKmer<SIZE> kmer;
		std::unique_ptr<WritingGathererBin<out_kmcdb_value_type>> outGahtererBin = gatherer.getBin(taskData.binId);
		for (auto kmer_idx = binsIndicesForCorrection[taskData.binId]; bin->NextKmer(kmer, inMatrixEntry.data()); ++kmer_idx)
		{
			kmer.to_string(kmer_len, kmerSequence.data());

			if (statisticsToGeneration.normalize)
			{
				normalization.norm_entry(params.statisticsParams.normalizationMethod, inMatrixEntry, outNormEntry);

				if (umap)
				{
					for (size_t sample_id = 0; sample_id < outNormEntry.size(); ++sample_id)
						umap->add(sample_id, kmer_idx, outNormEntry[sample_id]);
				}

				normOutputBuffer->StoreKmer(kmerSequence, outNormEntry, StoreMethods::AsMatrixRow);
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
				const double entropy = entropyObj.entropy_n(
					outNormEntry.begin(),
					num_samples);

				outStatsEntry[outStatsEntryIdx++] = entropy;
			}
			if (statisticsToGeneration.differentialAnalysis)
			{
				if (statisticsToGeneration.tTest)
				{
					const double tTestPValue = statistics.t_test_n(
						outNormEntry.begin(),
						differentialAnalysisPhenotype.begin(),
						num_samples).p_value;

					outStatsEntry[outStatsEntryIdx++] = tTestPValue;
				}
				if (statisticsToGeneration.snr)
				{
					const double snr = statistics.SNR_test_n(
						outNormEntry.begin(),
						differentialAnalysisPhenotype.begin(),
						num_samples);

					outStatsEntry[outStatsEntryIdx++] = snr;
				}
				if (statisticsToGeneration.wilcoxonRankSum)
				{
					const double wilcoxonRankSumPValue = statistics.mann_whitney_U_test_n(
						outNormEntry.begin(),
						differentialAnalysisPhenotype.begin(),
						num_samples).p_value;

					outStatsEntry[outStatsEntryIdx++] = wilcoxonRankSumPValue;
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
					const double anovaPValue = scorer.anova_n(
						outNormEntry.begin(),
						differentialAnalysisPhenotype.begin(),
						differentialAnalysisNClasses,
						num_samples).p_value;

					outStatsEntry[outStatsEntryIdx++] = anovaPValue;
				}
			}
			assert(outStatsEntryIdx == statisticsToGeneration.nStatistics);

			++progress_bar_updater;

			outGahtererBin->writeKmer(outStatsEntry, kmer, kmerSequence, inMatrixEntry, &keepNLargestCollection);
		}
	}
	keepNLargestCollectionGlobal.Add(keepNLargestCollection);
}



template<unsigned SIZE>
void StatisticsGenerator::processEntriesWhenCorrection(KeepNLargestCollectionGlobal<SIZE>& keepNLargestCollectionGlobal,
	refresh::umap_direct<double>* umap)
{
	assert(statisticsToGeneration.differentialAnalysis);
	assert(params.statisticsParams.generateNormalization);

	KeepNLargestCollection<SIZE> keepNLargestCollection;

	initKeepNLargest(keepNLargestCollection);

	TaskData taskData;
	while (tasksPool.getTask(taskData))
	{
		auto bin = matrixReader->GetBin(taskData.binId);

		uint64_t outPValuesToCorrectIdx = binsIndicesForCorrection[taskData.binId];
		const uint64_t outPValuesToCorrectIdxEnd = binsIndicesForCorrection[taskData.binId + 1];

		std::unique_ptr<OutputBuffer> normOutputBuffer = std::make_unique<OutputBuffer>(*normWriter, getMaxNormLineLength());

		refresh::normalization_work<uint64_t, double> normalization;
		normalization.register_method(params.statisticsParams.normalizationMethod);
		normalization.set_no_series(params.mkmcParams.samples.size());
		normalization.deserialize(params.statisticsParams.normalizationMethod, normalizationData);

		normalization.initialize();

		refresh::correlation correlation;
		refresh::statistics_entropy entropyObj;
		refresh::statistical_test statistics;
		refresh::scorers scorer;

		std::vector<uint64_t> inMatrixEntry;
		std::vector<out_kmcdb_value_type> outNormEntry; // normalized counts
		std::vector<out_kmcdb_value_type> outStatsEntry; // statistics
		std::ptrdiff_t num_samples = static_cast<std::ptrdiff_t>(params.mkmcParams.samples.size());

		inMatrixEntry.resize(num_samples);
		outNormEntry.resize(num_samples);
		outStatsEntry.resize(statisticsToGeneration.nStatistics - statisticsToGeneration.nStatisticsWithPValues);

		ProgressBarUpdater progress_bar_updater(*progress_bar, (std::max)(1ull, progress_bar->GetTotal() / 100ull));

		auto kmer_len = params.stage1Params.GetKmerLen();
		std::string kmerSequence(kmer_len, ' ');

		kmcdb::CKmer<SIZE> kmer;
		std::unique_ptr<WritingGathererBin<out_kmcdb_value_type>> outGathererBin = gatherer.getBin(taskData.binId, false, true);
		while (bin->NextKmer(kmer, inMatrixEntry.data()))
		{
			kmer.to_string(kmer_len, kmerSequence.data());

			normalization.norm_entry(params.statisticsParams.normalizationMethod, inMatrixEntry, outNormEntry);

			if (umap)
			{
				for (size_t sample_id = 0; sample_id < outNormEntry.size(); ++sample_id)
					umap->add(sample_id, outPValuesToCorrectIdx, outNormEntry[sample_id]);
			}

			normOutputBuffer->StoreKmer(kmerSequence, outNormEntry, StoreMethods::AsMatrixRow);

			size_t outStatsIdx = 0;
			size_t outPValuesToCorrectAlg = 0;

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
				const double entropy = entropyObj.entropy_n(
					outNormEntry.begin(),
					num_samples);

				outStatsEntry[outStatsIdx++] = entropy;
			}

			if (statisticsToGeneration.tTest)
			{
				const double tTestPValue = statistics.t_test_n(
					outNormEntry.begin(),
					differentialAnalysisPhenotype.begin(),
					num_samples).p_value;

				pValuesData[outPValuesToCorrectAlg++][outPValuesToCorrectIdx] = tTestPValue;
			}
			if (statisticsToGeneration.snr)
			{
				const double snr = statistics.SNR_test_n(
					outNormEntry.begin(),
					differentialAnalysisPhenotype.begin(),
					num_samples);

				outStatsEntry[outStatsIdx++] = snr;
			}
			if (statisticsToGeneration.wilcoxonRankSum)
			{
				const double wilcoxonRankSumPValue = statistics.mann_whitney_U_test_n(
					outNormEntry.begin(),
					differentialAnalysisPhenotype.begin(),
					num_samples).p_value;

				pValuesData[outPValuesToCorrectAlg++][outPValuesToCorrectIdx] = wilcoxonRankSumPValue;
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
				const double anovaPValue = scorer.anova_n(
					outNormEntry.begin(),
					differentialAnalysisPhenotype.begin(),
					differentialAnalysisNClasses,
					num_samples).p_value;

				pValuesData[outPValuesToCorrectAlg++][outPValuesToCorrectIdx] = anovaPValue;
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
	std::ptrdiff_t num_samples = static_cast<std::ptrdiff_t>(params.mkmcParams.samples.size());

	std::vector<uint64_t> inMatrixEntry;
	inMatrixEntry.resize(num_samples);

	TaskData taskData;
	while (tasksPool.getTask(taskData))
	{
		auto bin = matrixReader->GetBin(taskData.binId);

		uint64_t outPValuesToCorrectIdx = binsIndicesForCorrection[taskData.binId];
		const uint64_t outPValuesToCorrectIdxEnd = binsIndicesForCorrection[taskData.binId + 1];

		std::vector<out_kmcdb_value_type> outStatsEntry; // statistics

		outStatsEntry.resize(statisticsToGeneration.nStatistics);

		ProgressBarUpdater progress_bar_updater(*progress_bar, (std::max)(1ull, progress_bar->GetTotal() / 100ull));

		auto kmer_len = params.stage1Params.GetKmerLen();
		std::string kmerSequence(kmer_len, ' ');

		kmcdb::CKmer<SIZE> kmer;
		std::unique_ptr<WritingGathererBin<out_kmcdb_value_type>> outGathererBin = gatherer.getBin(taskData.binId, true, false);
		while (bin->NextKmer(kmer, inMatrixEntry.data()))
		{
			kmer.to_string(kmer_len, kmerSequence.data());

			size_t outStatsIdx = 0;
			if (statisticsToGeneration.differentialAnalysis)
			{
				if (statisticsToGeneration.tTest)
				{
					outStatsEntry[outStatsIdx] = pValuesCorrectedData[outStatsIdx][outPValuesToCorrectIdx];
					++outStatsIdx;
				}
				if (statisticsToGeneration.wilcoxonRankSum)
				{
					outStatsEntry[outStatsIdx] = pValuesCorrectedData[outStatsIdx][outPValuesToCorrectIdx];
					++outStatsIdx;
				}
				if (statisticsToGeneration.anova)
				{
					outStatsEntry[outStatsIdx] = pValuesCorrectedData[outStatsIdx][outPValuesToCorrectIdx];
					++outStatsIdx;
				}
			}
			assert(outStatsIdx == statisticsToGeneration.nStatisticsWithPValues);

			++outPValuesToCorrectIdx;
			++progress_bar_updater;

			outGathererBin->writeKmer(outStatsEntry, kmer, kmerSequence, inMatrixEntry);
		}
		assert(outPValuesToCorrectIdx == outPValuesToCorrectIdxEnd);
	}
}
