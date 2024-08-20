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

	const std::vector<int64_t>& correlationPhenotype;
	const std::vector<uint32_t>& differentialAnalysisPhenotype;
	size_t differentialAnalysisNClasses;

	std::vector<std::vector<out_kmcdb_value_type>> pValuesData; // first index: algorithm, second: entries
	std::vector<std::vector<out_kmcdb_value_type>> pValuesCorrectedData; // first index: algorithm, second: entries

	//mkokot_TODO: I am using this also for umap, so maybe change a name of this variable, like binsOffsets
	std::vector<uint64_t> binsIndicesForCorrection; // (of size no. of bins + 1) contains indices of first entries for every bin; the last element contains number of all the entries

	StatisticsToGeneration statisticsToGeneration;

	size_t getMaxNormLineLength() const
	{
		//@Maciej: I have added 1 after k-mer len as a separator
		return params.stage1Params.GetKmerLen() + 1 + params.mkmcParams.samples.size() * (refresh::numeric_conversion_max_length<out_kmcdb_value_type>() + 1);
	}

	template<unsigned SIZE>
	void initKeepNLargest(KeepNLargestCollection<SIZE>& keepNLargestCollection);

	template<unsigned SIZE>
	void processEntries(KeepNLargestCollectionGlobal<SIZE>& keepNLargestCollectionGlobal,
		refresh::umap_direct<double>* umap);

	void gatherPValuesEntriesToCorrection();

	void correctPValuesEntries();

	template<unsigned SIZE>
	void processEntriesAfterCorrection(KeepNLargestCollectionGlobal<SIZE>& keepNLargestCollectionGlobal,
		refresh::umap_direct<double>* umap);

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
		std::vector<out_kmcdb_value_type> outEntry; // normalized counts and statistics
		std::ptrdiff_t num_samples = static_cast<std::ptrdiff_t>(params.mkmcParams.samples.size());

		inMatrixEntry.resize(num_samples);
		outEntry.resize(statisticsToGeneration.nResults);

		ProgressBarUpdater progress_bar_updater(*progress_bar, (std::max)(1ull, progress_bar->GetTotal() / 100ull));

		auto kmer_len = params.stage1Params.GetKmerLen();
		std::string kmerSequence(kmer_len, ' ');

		kmcdb::CKmer<SIZE> kmer;
		std::unique_ptr<WritingGathererBin<out_kmcdb_value_type>> outBin = gatherer.getBin(taskData.binId);
		for (auto kmer_idx = binsIndicesForCorrection[taskData.binId]; bin->NextKmer(kmer, inMatrixEntry.data()); ++kmer_idx)
		{
			kmer.to_string(kmer_len, kmerSequence.data());

			if (statisticsToGeneration.normalize)
			{
				normalization.norm_entry(params.statisticsParams.normalizationMethod, inMatrixEntry, outEntry);

				if (umap)
				{
					for (size_t sample_id = 0; sample_id < outEntry.size(); ++sample_id)
						umap->add(sample_id, kmer_idx, outEntry[sample_id]);
				}

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

			outBin->writeKmer(outEntry, kmer, kmerSequence, inMatrixEntry, keepNLargestCollection);
		}
	}
	keepNLargestCollectionGlobal.Add(keepNLargestCollection);
}

template<unsigned SIZE>
void StatisticsGenerator::processEntriesAfterCorrection(KeepNLargestCollectionGlobal<SIZE>& keepNLargestCollectionGlobal,
	refresh::umap_direct<double>* umap)
{
	KeepNLargestCollection<SIZE> keepNLargestCollection;

	initKeepNLargest(keepNLargestCollection);

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

		kmcdb::CKmer<SIZE> kmer;
		std::unique_ptr<WritingGathererBin<out_kmcdb_value_type>> outBin = gatherer.getBin(taskData.binId);
		for (auto kmer_idx = binsIndicesForCorrection[taskData.binId]; bin->NextKmer(kmer, inMatrixEntry.data()); ++kmer_idx)
		{
			kmer.to_string(kmer_len, kmerSequence.data());

			if (statisticsToGeneration.normalize)
			{
				normalization.norm_entry(params.statisticsParams.normalizationMethod, inMatrixEntry, outEntry);

				if (umap)
				{
					for (size_t sample_id = 0; sample_id < outEntry.size(); ++sample_id)
						umap->add(sample_id, kmer_idx, outEntry[sample_id]);
				}

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

			outBin->writeKmer(outEntry, kmer, kmerSequence, inMatrixEntry, keepNLargestCollection);
		}

		assert(dataIdx == dataIdxEnd);
	}
	keepNLargestCollectionGlobal.Add(keepNLargestCollection);
}
