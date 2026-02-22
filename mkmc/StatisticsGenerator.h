#pragma once

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



template<unsigned SIZE, typename out_kmcdb_value_type, typename cnt_value_type>
class CVGenerator
{
	const bool cv;
	const std::vector<out_kmcdb_value_type>& wholeCorrelationPhenotype;

	const size_t nSamples;
	const std::vector<size_t> samplesToBeTestOrder;

	const size_t nTestSamples;
	const size_t nFolds;
	const size_t nTrainingSamples;

	std::vector<out_kmcdb_value_type> entry;
	std::vector<out_kmcdb_value_type> correlationPhenotype; // Actually we can pregenerate all the possible sequences and store them, but updating this vector for every iteration shouldn't be a huge cost

	// If not empty, compute the correlation
	std::vector<out_kmcdb_value_type> outPearsonStatsEntry, outSpearmanStatsEntry, outKendallStatsEntry;

	KeepNLargestCollectionCV<SIZE, out_kmcdb_value_type, cnt_value_type> keepNLargestCollection;
	KeepNLargestCollectionGlobal<SIZE, out_kmcdb_value_type, cnt_value_type, KeepNLargestCollectionCV<SIZE, out_kmcdb_value_type, cnt_value_type>> &keepNLargestCollectionCVGlobal;

	refresh::correlation& correlation;

public:
	CVGenerator(Params& params,
		const std::vector<out_kmcdb_value_type>& wholeCorrelationPhenotype,
		size_t numSamples,
		bool pearson,
		bool spearman,
		bool kendall,
		KeepNLargestCollectionGlobal<SIZE, out_kmcdb_value_type, cnt_value_type, KeepNLargestCollectionCV<SIZE, out_kmcdb_value_type, cnt_value_type>>&keepNLargestCollectionCVGlobal,
		refresh::correlation& correlation) :
		cv(params.statisticsParams.cvParams.cv),
		wholeCorrelationPhenotype(wholeCorrelationPhenotype),
		nSamples(numSamples),
		samplesToBeTestOrder(params.statisticsParams.cvParams.samplesToBeTestOrder),
		nTestSamples(params.statisticsParams.cvParams.nTestSamples),
		nFolds(numSamples / nTestSamples),
		nTrainingSamples(numSamples - nTestSamples),
		entry(nTrainingSamples),
		correlationPhenotype(nTrainingSamples),
		keepNLargestCollection(samplesToBeTestOrder,
			params.statisticsParams.nTop,
			nFolds,
			pearson,
			spearman,
			kendall),
		keepNLargestCollectionCVGlobal(keepNLargestCollectionCVGlobal),
		correlation(correlation)
	{
		if (params.statisticsParams.cvParams.cv)
		{
			assert(pearson || spearman || kendall);
			assert(params.statisticsParams.generateNormalization);

			if (pearson)
				outPearsonStatsEntry.resize(nFolds);
			if (spearman)
				outSpearmanStatsEntry.resize(nFolds);
			if (kendall)
				outKendallStatsEntry.resize(nFolds);
		}
	}

	~CVGenerator() {
		keepNLargestCollectionCVGlobal.Add(keepNLargestCollection);
	}

	void correlationCV(const kmcdb::CKmer<SIZE>& kmer, const std::string& kmerSeq, const std::vector<cnt_value_type>& wholeInputMatrixEntry, const std::vector<out_kmcdb_value_type>& wholeEntry);
};



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
	void determineBinsSizes();
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
	void processEntries(KeepNLargestCollectionGlobal<SIZE, out_kmcdb_value_type, cnt_value_type, KeepNLargestCollection<SIZE, out_kmcdb_value_type, cnt_value_type>>& keepNLargestCollectionGlobal,
		KeepNLargestCollectionGlobal<SIZE, out_kmcdb_value_type, cnt_value_type, KeepNLargestCollectionCV<SIZE, out_kmcdb_value_type, cnt_value_type>>& keepNLargestCollectionCVGlobal,
		DimensionalityReduction& dimensionalityReduction);

	template<unsigned SIZE>
	void processEntriesWhenCorrection(KeepNLargestCollectionGlobal<SIZE, out_kmcdb_value_type, cnt_value_type, KeepNLargestCollection<SIZE, out_kmcdb_value_type, cnt_value_type>>& keepNLargestCollectionGlobal,
		KeepNLargestCollectionGlobal<SIZE, out_kmcdb_value_type, cnt_value_type, KeepNLargestCollectionCV<SIZE, out_kmcdb_value_type, cnt_value_type>>& keepNLargestCollectionCVGlobal,
		DimensionalityReduction& dimensionalityReduction);

	void correctPValuesEntries();

	template<unsigned SIZE>
	void safeCorrectedPValuesEntries();

	bool readNormalizationData();

public:
	StatisticsGenerator(Params& params);

	void generateStatisticsParallel();
};


template<unsigned SIZE, typename out_kmcdb_value_type, typename cnt_value_type>
void CVGenerator<SIZE, out_kmcdb_value_type, cnt_value_type>::correlationCV(const kmcdb::CKmer<SIZE>& kmer, const std::string& kmerSeq, const std::vector<cnt_value_type>& wholeInputMatrixEntry, const std::vector<out_kmcdb_value_type>& wholeEntry)
{
	if (!cv)
		return;

	// For wholeEntry = ABCDEFGH, nTestSamples = 2, and samplesToBeTestOrder = 01234567
	// treat as test samples 0 and 1, then 2 and 3...
	// Entry will contain training subsequences of wholeEntry samples counts of size nTestSamples:
	// CDEFGH
	// ABEFGH
	// ABCDGH
	// ABCDEF
	// For wholeEntry = ABCDEFGH, nTestSamples = 2, and samplesToBeTestOrder = 57041326:
	// AEBDCG
	// FHBDCG
	// FHAECG
	// FHAEBD
	// It is required that samplesToBeTestOrder contains permutation of values 0...samples-1, but providing that
	// every subsequence 0...nTestSamples-1, nTestSamples...2*nTestSamples-1, ..., nSamples-nTestSamples...nSamples-1 is sorted.
	// E.g. 57041326 is OK, but 75041326 is wrong for nTestSamples.
	
	size_t outStatsEntryIdx = 0;
	size_t iTestSamples = 0;
	for (size_t iFold = 0; iFold < nFolds; ++iFold)
	{
		size_t iCurrentFoldTestSamples = 0;
		for (size_t iSamples = 0; iSamples < nSamples; ++iSamples)
		{
			// Include to training set all samples, except the ones present in a proper subsequence of samplesToBeTestOrder
			assert(iTestSamples < nTestSamples * (iFold + 1) == iCurrentFoldTestSamples < nTestSamples);
			if (iCurrentFoldTestSamples < nTestSamples && iSamples == samplesToBeTestOrder[iTestSamples])
			{
				++iTestSamples;
				++iCurrentFoldTestSamples;
				continue;
			}
			entry[iSamples - iCurrentFoldTestSamples] = wholeEntry[iSamples];
			correlationPhenotype[iSamples - iCurrentFoldTestSamples] = wholeCorrelationPhenotype[iSamples];
		}

		if (!outPearsonStatsEntry.empty())
		{
			const double pearson = refresh::correlation::pearson_n(
				entry.begin(),
				correlationPhenotype.begin(),
				nTrainingSamples);

			outPearsonStatsEntry[outStatsEntryIdx] = pearson;
		}
		if (!outSpearmanStatsEntry.empty())
		{
			const double spearman = correlation.spearman_n(
				entry.begin(),
				correlationPhenotype.begin(),
				nTrainingSamples);

			outSpearmanStatsEntry[outStatsEntryIdx] = spearman;
		}
		if (!outKendallStatsEntry.empty())
		{
			const double kendall = refresh::correlation::kendall_tau_n(
				entry.begin(),
				correlationPhenotype.begin(),
				nTrainingSamples);

			outKendallStatsEntry[outStatsEntryIdx] = kendall;
		}

		++outStatsEntryIdx;
	}

	keepNLargestCollection.addPearson(kmerSeq, kmer, outPearsonStatsEntry, wholeInputMatrixEntry);
	keepNLargestCollection.addSpearman(kmerSeq, kmer, outSpearmanStatsEntry, wholeInputMatrixEntry);
	keepNLargestCollection.addKendall(kmerSeq, kmer, outKendallStatsEntry, wholeInputMatrixEntry);
}



template<unsigned SIZE>
void StatisticsGenerator::processEntries(KeepNLargestCollectionGlobal<SIZE, out_kmcdb_value_type, cnt_value_type, KeepNLargestCollection<SIZE, out_kmcdb_value_type, cnt_value_type>>& keepNLargestCollectionGlobal,
	KeepNLargestCollectionGlobal<SIZE, out_kmcdb_value_type, cnt_value_type, KeepNLargestCollectionCV<SIZE, out_kmcdb_value_type, cnt_value_type>>& keepNLargestCollectionCVGlobal,
	DimensionalityReduction& dimensionalityReduction)
{
	const std::size_t num_samples = params.mkmcParams.samples.size();

	KeepNLargestCollection<SIZE, out_kmcdb_value_type, cnt_value_type> keepNLargestCollection(params.statisticsParams.nTop,
		statisticsToGeneration.pearson,
		statisticsToGeneration.spearman,
		statisticsToGeneration.kendall,
		statisticsToGeneration.entropy,
		statisticsToGeneration.snr,
		statisticsToGeneration.unnormalizedSnr,
		statisticsToGeneration.dids);

	refresh::correlation correlation;
	refresh::statistics_entropy entropyObj;
	refresh::statistical_test statistics;
	refresh::scorers scorer;

	CVGenerator<SIZE, out_kmcdb_value_type, cnt_value_type> cvGenerator(params,
		correlationPhenotype,
		num_samples,
		statisticsToGeneration.pearson,
		statisticsToGeneration.spearman,
		statisticsToGeneration.kendall,
		keepNLargestCollectionCVGlobal,
		correlation);

	std::vector<cnt_value_type> inMatrixEntry;
	std::vector<out_kmcdb_value_type> outNormEntry; // normalized counts
	std::vector<out_kmcdb_value_type> outStatsEntry; // statistics
	inMatrixEntry.resize(num_samples);
	outNormEntry.resize(num_samples);
	outStatsEntry.resize(statisticsToGeneration.nStatistics + statisticsToGeneration.nAdditionalValuesOfCorrectedStats);

	const auto kmer_len = params.stage1Params.GetKmerLen();
	std::string kmerSequence(kmer_len, ' ');

	TaskData taskData;
	while (tasksPool.getTask(taskData))
	{
		auto bin = matrixReader->GetBin(taskData.binId);

		std::unique_ptr<MatrixOutputBuffer<out_kmcdb_value_type>> normOutputBuffer;

		refresh::normalization_work<cnt_value_type, out_kmcdb_value_type> normalization;
		if (statisticsToGeneration.normalize)
		{
			normalization.register_method(params.statisticsParams.normalizationMethod);
			normalization.set_no_series(params.mkmcParams.samples.size());
			normalization.deserialize(params.statisticsParams.normalizationMethod, normalizationData);

			normalization.initialize();

			if (statisticsToGeneration.saveNormalization)
				normOutputBuffer = std::make_unique<MatrixOutputBuffer<out_kmcdb_value_type>>(*normWriter, kmer_len, num_samples);
		}

		ProgressBarUpdater progress_bar_updater(*progress_bar, (std::max)(1ull, progress_bar->GetTotal() / 100ull));

		kmcdb::CKmer<SIZE> kmer;
		std::unique_ptr<WritingGathererBin<out_kmcdb_value_type, cnt_value_type>> outGahtererBin = gatherer.getBin(taskData.binId);

		std::vector<double> inMatrixEntryScaled(num_samples); // for t-test
		for (auto kmer_idx = binsOffsets[taskData.binId]; bin->NextKmer(kmer, inMatrixEntry.data()); ++kmer_idx)
		{
			kmer.to_string(kmer_len, kmerSequence.data());

			if (statisticsToGeneration.normalize)
			{
				normalization.norm_entry(params.statisticsParams.normalizationMethod, inMatrixEntry, outNormEntry);

				dimensionalityReduction.add(kmer_idx, outNormEntry);

				if (statisticsToGeneration.saveNormalization)
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
					typedef StatisticsParams::DIDSMode DIDSMode;
					double dids;
					switch (params.statisticsParams.didsMode)
					{
					case DIDSMode::sqrt:
						dids = scorer.dids_n(
							outNormEntry.begin(),
							differentialAnalysisPhenotype.begin(),
							differentialAnalysisNClasses,
							num_samples);
						break;
					case DIDSMode::quadratic:
						dids = scorer.dids_quadratic_n(
							outNormEntry.begin(),
							differentialAnalysisPhenotype.begin(),
							differentialAnalysisNClasses,
							num_samples);
						break;
					case DIDSMode::tanh:
						dids = scorer.dids_tanh_n(
							outNormEntry.begin(),
							differentialAnalysisPhenotype.begin(),
							differentialAnalysisNClasses,
							num_samples);
						break;
					default:
						assert(false);
					}

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

			cvGenerator.correlationCV(kmer, kmerSequence, inMatrixEntry, outNormEntry);

			outGahtererBin->writeKmer(outStatsEntry, kmer, kmerSequence, inMatrixEntry, taskData.binId, kmer_idx - binsOffsets[taskData.binId], &keepNLargestCollection);
			++progress_bar_updater;
		}
	}
	keepNLargestCollectionGlobal.Add(keepNLargestCollection);
}



template<unsigned SIZE>
void StatisticsGenerator::processEntriesWhenCorrection(KeepNLargestCollectionGlobal<SIZE, out_kmcdb_value_type, cnt_value_type, KeepNLargestCollection<SIZE, out_kmcdb_value_type, cnt_value_type>>& keepNLargestCollectionGlobal,
	KeepNLargestCollectionGlobal<SIZE, out_kmcdb_value_type, cnt_value_type, KeepNLargestCollectionCV<SIZE, out_kmcdb_value_type, cnt_value_type>>& keepNLargestCollectionCVGlobal,
	DimensionalityReduction& dimensionalityReduction)
{
	assert(statisticsToGeneration.differentialAnalysis);

	const std::size_t num_samples = params.mkmcParams.samples.size();

	KeepNLargestCollection<SIZE, out_kmcdb_value_type, cnt_value_type> keepNLargestCollection(params.statisticsParams.nTop,
		statisticsToGeneration.pearson,
		statisticsToGeneration.spearman,
		statisticsToGeneration.kendall,
		statisticsToGeneration.entropy,
		statisticsToGeneration.snr,
		statisticsToGeneration.unnormalizedSnr,
		statisticsToGeneration.dids);

	refresh::correlation correlation;
	refresh::statistics_entropy entropyObj;
	refresh::statistical_test statistics;
	refresh::scorers scorer;

	CVGenerator<SIZE, out_kmcdb_value_type, cnt_value_type> cvGenerator(params,
		correlationPhenotype,
		num_samples,
		statisticsToGeneration.pearson,
		statisticsToGeneration.spearman,
		statisticsToGeneration.kendall,
		keepNLargestCollectionCVGlobal,
		correlation);

	std::vector<cnt_value_type> inMatrixEntry;
	std::vector<out_kmcdb_value_type> outNormEntry; // normalized counts
	std::vector<out_kmcdb_value_type> outStatsEntry; // statistics
	inMatrixEntry.resize(num_samples);
	outNormEntry.resize(num_samples);
	outStatsEntry.resize(statisticsToGeneration.nStatistics - statisticsToGeneration.nStatisticsWithPValues);

	const auto kmer_len = params.stage1Params.GetKmerLen();
	std::string kmerSequence(kmer_len, ' ');

	TaskData taskData;
	while (tasksPool.getTask(taskData))
	{
		auto bin = matrixReader->GetBin(taskData.binId);

		uint64_t outPValuesToCorrectIdx = binsOffsets[taskData.binId];
		const uint64_t outPValuesToCorrectIdxEnd = binsOffsets[taskData.binId + 1];

		std::unique_ptr<MatrixOutputBuffer<out_kmcdb_value_type>> normOutputBuffer;

		refresh::normalization_work<cnt_value_type, out_kmcdb_value_type> normalization;
		if (statisticsToGeneration.normalize)
		{
			normalization.register_method(params.statisticsParams.normalizationMethod);
			normalization.set_no_series(params.mkmcParams.samples.size());
			normalization.deserialize(params.statisticsParams.normalizationMethod, normalizationData);

			normalization.initialize();

			if (statisticsToGeneration.saveNormalization)
				normOutputBuffer = std::make_unique<MatrixOutputBuffer<out_kmcdb_value_type>>(*normWriter, kmer_len, num_samples);
		}

		ProgressBarUpdater progress_bar_updater(*progress_bar, (std::max)(1ull, progress_bar->GetTotal() / 100ull));

		kmcdb::CKmer<SIZE> kmer;
		std::unique_ptr<WritingGathererBin<out_kmcdb_value_type, cnt_value_type>> outGathererBin = gatherer.getBin(taskData.binId, false, true);

		std::vector<double> inMatrixEntryScaled(num_samples); // for t-test
		while (bin->NextKmer(kmer, inMatrixEntry.data()))
		{
			kmer.to_string(kmer_len, kmerSequence.data());

			if (statisticsToGeneration.normalize)
			{
				normalization.norm_entry(params.statisticsParams.normalizationMethod, inMatrixEntry, outNormEntry);

				dimensionalityReduction.add(outPValuesToCorrectIdx, outNormEntry);

				if (statisticsToGeneration.saveNormalization)
					normOutputBuffer->StoreKmer(kmerSequence, outNormEntry);
			}

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
				typedef StatisticsParams::DIDSMode DIDSMode;
				double dids;
				switch (params.statisticsParams.didsMode)
				{
				case DIDSMode::sqrt:
					dids = scorer.dids_n(
						outNormEntry.begin(),
						differentialAnalysisPhenotype.begin(),
						differentialAnalysisNClasses,
						num_samples);
					break;
				case DIDSMode::quadratic:
					dids = scorer.dids_quadratic_n(
						outNormEntry.begin(),
						differentialAnalysisPhenotype.begin(),
						differentialAnalysisNClasses,
						num_samples);
					break;
				case DIDSMode::tanh:
					dids = scorer.dids_tanh_n(
						outNormEntry.begin(),
						differentialAnalysisPhenotype.begin(),
						differentialAnalysisNClasses,
						num_samples);
					break;
				default:
					assert(false);
				}

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

			cvGenerator.correlationCV(kmer, kmerSequence, inMatrixEntry, outNormEntry);

			outGathererBin->writeKmer(outStatsEntry, kmer, kmerSequence, inMatrixEntry, taskData.binId, outPValuesToCorrectIdx - binsOffsets[taskData.binId], &keepNLargestCollection);
			++outPValuesToCorrectIdx;
			++progress_bar_updater;
		}

		assert(outPValuesToCorrectIdx == outPValuesToCorrectIdxEnd);
	}
	keepNLargestCollectionGlobal.Add(keepNLargestCollection);
}



template<unsigned SIZE>
void StatisticsGenerator::safeCorrectedPValuesEntries()
{
	assert(statisticsToGeneration.differentialAnalysis); // Currently, nTestSamples-values are generated for DA statistics only
	const std::size_t num_samples = params.mkmcParams.samples.size();

	std::vector<cnt_value_type> inMatrixEntry;
	std::vector<out_kmcdb_value_type> outStatsEntry; // statistics
	inMatrixEntry.resize(num_samples);
	outStatsEntry.resize(statisticsToGeneration.nStatistics + statisticsToGeneration.nAdditionalValuesOfCorrectedStats);

	const auto kmer_len = params.stage1Params.GetKmerLen();
	std::string kmerSequence(kmer_len, ' ');

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
			kmer.to_string(kmer_len, kmerSequence.data());

			size_t outPvaluesIdx = 0;
			size_t outAdditionalValuesIdx = 0;
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
			assert(outPvaluesIdx + outAdditionalValuesIdx == statisticsToGeneration.nStatisticsWithPValues + statisticsToGeneration.nAdditionalValuesOfCorrectedStats);

			outGathererBin->writeKmer(outStatsEntry, kmer, kmerSequence, inMatrixEntry, taskData.binId, outPValuesToCorrectIdx - binsOffsets[taskData.binId]);
			++outPValuesToCorrectIdx;
			++progress_bar_updater;
		}
		assert(outPValuesToCorrectIdx == outPValuesToCorrectIdxEnd);
	}
}
