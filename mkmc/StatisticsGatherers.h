#pragma once

#include <memory>
#include <vector>
#include <string>

#include "Merger.h"
#include "parameters.h"
#include "TextFileWritingUtilities.h"
#include "kmcdb/kmcdb.h"
#include "KeepNLargests.h"
#include "refresh/statistics/lib/statistics_umap.h"


struct StatisticsToGeneration
{
	bool normalize = false;
	bool saveNormalization = false;

	bool pearson = false;
	bool spearman = false;
	bool kendall = false;

	bool entropy = false;

	bool differentialAnalysis = false; // logical sum of the following ones

	bool tTest = false;
	bool snr = false;
	bool unnormalizedSnr = false;
	bool wrs = false;

	bool dids = false;
	bool anova = false;

	uint32_t nStatistics = 0;
	uint32_t nStatisticsWithPValues = 0;
	uint32_t nAdditionalValuesOfCorrectedStats = 0; // T-Test = 2, WRS = 2, ANOVA = 1
};



template<typename Statistics_T, typename VALUE_T>
class WritingGatherer;



template<typename Statistics_T, typename VALUE_T>
class WritingGathererBin
{
	friend WritingGatherer<Statistics_T, VALUE_T>;
	WritingGatherer<Statistics_T, VALUE_T>& mainWritingGatherer;

	const uint32_t kmerLength;
	const uint32_t numSamples;
	const double maxCorrectedPval;

	struct {
		std::unique_ptr<MatrixOutputBuffer<Statistics_T>> pearson;
		std::unique_ptr<MatrixOutputBuffer<Statistics_T>> spearman;
		std::unique_ptr<MatrixOutputBuffer<Statistics_T>> kendall;

		std::unique_ptr<MatrixOutputBuffer<Statistics_T>> entropy;

		std::unique_ptr<MatrixOutputBuffer<Statistics_T>> tTest;
		std::unique_ptr<MatrixOutputBuffer<Statistics_T>> tTestSignificant;
		std::unique_ptr<MatrixOutputBuffer<VALUE_T>> tTestSignificantCntMatrix;
		std::unique_ptr<FastaOutputBuffer> tTestSignificantFasta;

		std::unique_ptr<MatrixOutputBuffer<Statistics_T>> snr;
		std::unique_ptr<MatrixOutputBuffer<Statistics_T>> unnormalizedSnr;

		std::unique_ptr<MatrixOutputBuffer<Statistics_T>> wrs;
		std::unique_ptr<MatrixOutputBuffer<Statistics_T>> wrsSignificant;
		std::unique_ptr<MatrixOutputBuffer<VALUE_T>> wrsSignificantCntMatrix;
		std::unique_ptr<FastaOutputBuffer> wrsSignificantFasta;

		std::unique_ptr<MatrixOutputBuffer<Statistics_T>> dids;

		std::unique_ptr<MatrixOutputBuffer<Statistics_T>> anova;
		std::unique_ptr<MatrixOutputBuffer<Statistics_T>> anovaSignificant;
		std::unique_ptr<MatrixOutputBuffer<VALUE_T>> anovaSignificantCntMatrix;
		std::unique_ptr<FastaOutputBuffer> anovaSignificantFasta;
	} outputBuffers;

	WritingGathererBin(WritingGatherer<Statistics_T, VALUE_T>& mainWritingGatherer,
		const uint32_t kmerLength,
		uint32_t binId,
		uint32_t numSamples,
		double maxCorrectedPval,
		bool gatherCorrected,
		bool gatherNotCorrected);

public:
	template<unsigned SIZE>
	void writeKmer(
		const std::vector<Statistics_T>& outEntry, // outEntry - consecutive statistics values in a proper order
		const kmcdb::CKmer<SIZE>& kmer,
		const std::string& kmerSeq,
		const std::vector<VALUE_T>& original_counts,
		uint32_t binId,
		uint64_t kmerIdInBin,
		KeepNLargestCollection<SIZE, Statistics_T, VALUE_T>* keepNLargestCollection = nullptr);

};



template<typename Statistics_T, typename VALUE_T>
class WritingGatherer
{
	friend WritingGathererBin<Statistics_T, VALUE_T>;
	const Params& params;

	struct
	{
		std::unique_ptr<TextFileWriter> pearson;
		std::unique_ptr<TextFileWriter> spearman;
		std::unique_ptr<TextFileWriter> kendall;

		std::unique_ptr<TextFileWriter> entropy;

		std::unique_ptr<TextFileWriter> tTest;
		std::unique_ptr<TextFileWriter> tTestSignificant;
		std::unique_ptr<TextFileWriter> tTestSignificantCntMatrix;
		std::unique_ptr<TextFileWriter> tTestSignificantFasta;

		std::unique_ptr<TextFileWriter> snr;
		std::unique_ptr<TextFileWriter> unnormalizedSnr;

		std::unique_ptr<TextFileWriter> wrs;
		std::unique_ptr<TextFileWriter> wrsSignificant;
		std::unique_ptr<TextFileWriter> wrsSignificantCntMatrix;
		std::unique_ptr<TextFileWriter> wrsSignificantFasta;

		std::unique_ptr<TextFileWriter> dids;

		std::unique_ptr<TextFileWriter> anova;
		std::unique_ptr<TextFileWriter> anovaSignificant;
		std::unique_ptr<TextFileWriter> anovaSignificantCntMatrix;
		std::unique_ptr<TextFileWriter> anovaSignificantFasta;
	} writers;

	const StatisticsToGeneration& statisticsToGeneration;

public:
	WritingGatherer(const Params& params, const StatisticsToGeneration& statisticsToGeneration) :
		params(params),
		statisticsToGeneration(statisticsToGeneration)
	{}

	void initWriting(
		const std::unique_ptr<kmcdb::MetadataReader>& matrixMetadataReader,
		const std::vector<std::string>& cnt_matrix_output_header);

	std::unique_ptr<WritingGathererBin<Statistics_T, VALUE_T>> getBin(
		uint32_t binId,
		bool gatherCorrected = true,
		bool gatherNotCorrected = true)
	{
		// explicit new calling is necessary due to constructor's privacy
		return std::unique_ptr<WritingGathererBin<Statistics_T, VALUE_T>>(new WritingGathererBin<Statistics_T, VALUE_T>(
			*this,
			params.stage1Params.GetKmerLen(),
			binId,
			static_cast<uint32_t>(params.mkmcParams.samples.size()),
			params.statisticsParams.maxCorrectedPval,
			gatherCorrected,
			gatherNotCorrected
		));
	}
};



template<typename Statistics_T, typename VALUE_T>
WritingGathererBin<Statistics_T, VALUE_T>::WritingGathererBin(
	WritingGatherer<Statistics_T, VALUE_T>& mainWritingGatherer,
	const uint32_t kmerLength,
	uint32_t binId,
	uint32_t numSamples,
	double maxCorrectedPval,
	bool gatherCorrected,
	bool gatherNotCorrected) :
	mainWritingGatherer(mainWritingGatherer),
	kmerLength(kmerLength),
	numSamples(numSamples),
	maxCorrectedPval(maxCorrectedPval)
{
	assert((gatherCorrected && gatherNotCorrected) || mainWritingGatherer.params.statisticsParams.correctPvalues); // safe separately only if correction is performed
	assert(gatherCorrected || gatherNotCorrected);

	if (gatherNotCorrected && mainWritingGatherer.writers.pearson)
		outputBuffers.pearson = std::make_unique<MatrixOutputBuffer<Statistics_T>>(*mainWritingGatherer.writers.pearson, kmerLength, 1);
	if (gatherNotCorrected && mainWritingGatherer.writers.spearman)
		outputBuffers.spearman = std::make_unique<MatrixOutputBuffer<Statistics_T>>(*mainWritingGatherer.writers.spearman, kmerLength, 1);
	if (gatherNotCorrected && mainWritingGatherer.writers.kendall)
		outputBuffers.kendall = std::make_unique<MatrixOutputBuffer<Statistics_T>>(*mainWritingGatherer.writers.kendall, kmerLength, 1);

	if (gatherNotCorrected && mainWritingGatherer.writers.entropy)
		outputBuffers.entropy = std::make_unique<MatrixOutputBuffer<Statistics_T>>(*mainWritingGatherer.writers.entropy, kmerLength, 1);

	if (gatherCorrected && mainWritingGatherer.writers.tTest)
	{
		outputBuffers.tTest = std::make_unique<MatrixOutputBuffer<Statistics_T>>(*mainWritingGatherer.writers.tTest, kmerLength, 3);
		if (mainWritingGatherer.writers.tTestSignificant)
		{
			assert(mainWritingGatherer.writers.tTestSignificantCntMatrix);
			assert(mainWritingGatherer.writers.tTestSignificantFasta);

			outputBuffers.tTestSignificant = std::make_unique<MatrixOutputBuffer<Statistics_T>>(*mainWritingGatherer.writers.tTestSignificant, kmerLength, 3);
			outputBuffers.tTestSignificantCntMatrix = std::make_unique<MatrixOutputBuffer<VALUE_T>>(*mainWritingGatherer.writers.tTestSignificantCntMatrix, kmerLength, numSamples);
			outputBuffers.tTestSignificantFasta = std::make_unique<FastaOutputBuffer>(*mainWritingGatherer.writers.tTestSignificantFasta, kmerLength);
		}
	}
	if (gatherNotCorrected && mainWritingGatherer.writers.snr)
		outputBuffers.snr = std::make_unique<MatrixOutputBuffer<Statistics_T>>(*mainWritingGatherer.writers.snr, kmerLength, 1);
	if (gatherNotCorrected && mainWritingGatherer.writers.unnormalizedSnr)
		outputBuffers.unnormalizedSnr = std::make_unique<MatrixOutputBuffer<Statistics_T>>(*mainWritingGatherer.writers.unnormalizedSnr, kmerLength, 1);
	if (gatherCorrected && mainWritingGatherer.writers.wrs)
	{
		outputBuffers.wrs = std::make_unique<MatrixOutputBuffer<Statistics_T>>(*mainWritingGatherer.writers.wrs, kmerLength, 3);
		if (mainWritingGatherer.writers.wrsSignificant)
		{
			assert(mainWritingGatherer.writers.wrsSignificantCntMatrix);
			assert(mainWritingGatherer.writers.wrsSignificantFasta);

			outputBuffers.wrsSignificant = std::make_unique<MatrixOutputBuffer<Statistics_T>>(*mainWritingGatherer.writers.wrsSignificant, kmerLength, 3);
			outputBuffers.wrsSignificantCntMatrix = std::make_unique<MatrixOutputBuffer<VALUE_T>>(*mainWritingGatherer.writers.wrsSignificantCntMatrix, kmerLength, numSamples);
			outputBuffers.wrsSignificantFasta = std::make_unique<FastaOutputBuffer>(*mainWritingGatherer.writers.wrsSignificantFasta, kmerLength);
		}
	}

	if (gatherNotCorrected && mainWritingGatherer.writers.dids)
		outputBuffers.dids = std::make_unique<MatrixOutputBuffer<Statistics_T>>(*mainWritingGatherer.writers.dids, kmerLength, 1);
	if (gatherCorrected && mainWritingGatherer.writers.anova)
	{
		outputBuffers.anova = std::make_unique<MatrixOutputBuffer<Statistics_T>>(*mainWritingGatherer.writers.anova, kmerLength, 3);
		if (mainWritingGatherer.writers.anovaSignificant)
		{
			assert(mainWritingGatherer.writers.anovaSignificantCntMatrix);
			assert(mainWritingGatherer.writers.anovaSignificantFasta);

			outputBuffers.anovaSignificant = std::make_unique<MatrixOutputBuffer<Statistics_T>>(*mainWritingGatherer.writers.anovaSignificant, kmerLength, 3);
			outputBuffers.anovaSignificantCntMatrix = std::make_unique<MatrixOutputBuffer<VALUE_T>>(*mainWritingGatherer.writers.anovaSignificantCntMatrix, kmerLength, numSamples);
			outputBuffers.anovaSignificantFasta = std::make_unique<FastaOutputBuffer>(*mainWritingGatherer.writers.anovaSignificantFasta, kmerLength);
		}
	}
}



template<typename Statistics_T, typename VALUE_T>
template<unsigned SIZE>
void WritingGathererBin<Statistics_T, VALUE_T>::writeKmer(
	const std::vector<Statistics_T>& outEntry,
	const kmcdb::CKmer<SIZE>& kmer,
	const std::string& kmerSeq,
	const std::vector<VALUE_T>& original_counts,
	uint32_t binId,
	uint64_t kmerIdInBin,
	KeepNLargestCollection<SIZE, Statistics_T, VALUE_T>* keepNLargestCollection/* = nullptr*/)
{
	assert(keepNLargestCollection != nullptr || (!outputBuffers.pearson && !outputBuffers.spearman && !outputBuffers.kendall && !outputBuffers.entropy && !outputBuffers.snr && !outputBuffers.unnormalizedSnr && !outputBuffers.dids));
	size_t valuesIdx = 0;
	if (outputBuffers.pearson)
	{
		const auto value = outEntry[valuesIdx++];
		outputBuffers.pearson->StoreKmer(kmerSeq, value);
		keepNLargestCollection->addPearson(kmerSeq, kmer, value, original_counts);
	}
	if (outputBuffers.spearman)
	{
		const auto value = outEntry[valuesIdx++];
		outputBuffers.spearman->StoreKmer(kmerSeq, value);
		keepNLargestCollection->addSpearman(kmerSeq, kmer, value, original_counts);
	}
	if (outputBuffers.kendall)
	{
		const auto value = outEntry[valuesIdx++];
		outputBuffers.kendall->StoreKmer(kmerSeq, value);
		keepNLargestCollection->addKendall(kmerSeq, kmer, value, original_counts);
	}

	if (outputBuffers.entropy)
	{
		const auto value = outEntry[valuesIdx++];
		outputBuffers.entropy->StoreKmer(kmerSeq, value);
		keepNLargestCollection->addEntropy(kmerSeq, kmer, value, original_counts);
	}
	if (mainWritingGatherer.statisticsToGeneration.differentialAnalysis)
	{
		if (outputBuffers.tTest)
		{
			const auto pValue = outEntry[valuesIdx++];
			const auto df = outEntry[valuesIdx++];
			const auto statistic = outEntry[valuesIdx++];

			const std::vector<Statistics_T> valuesToSave{ pValue, df, statistic };

			outputBuffers.tTest->StoreKmer(kmerSeq, valuesToSave);

			if (outputBuffers.tTestSignificant)
			{
				if (pValue <= maxCorrectedPval)
				{
					outputBuffers.tTestSignificant->StoreKmer(kmerSeq, valuesToSave);
					outputBuffers.tTestSignificantCntMatrix->StoreKmer(kmerSeq, original_counts);
					outputBuffers.tTestSignificantFasta->StoreKmer(kmerSeq, binId, kmerIdInBin);
				}
			}
		}
		if (outputBuffers.snr)
		{
			const auto value = outEntry[valuesIdx++];
			outputBuffers.snr->StoreKmer(kmerSeq, value);
			keepNLargestCollection->addSnr(kmerSeq, kmer, value, original_counts);
		}
		if (outputBuffers.unnormalizedSnr)
		{
			const auto value = outEntry[valuesIdx++];
			outputBuffers.unnormalizedSnr->StoreKmer(kmerSeq, value);
			keepNLargestCollection->addUnnormalizedSnr(kmerSeq, kmer, value, original_counts);
		}
		if (outputBuffers.wrs)
		{
			const auto pValue = outEntry[valuesIdx++];
			const auto statisticU1 = outEntry[valuesIdx++];
			const auto statisticU2 = outEntry[valuesIdx++];

			const std::vector<Statistics_T> valuesToSave{ pValue, statisticU1, statisticU2 };

			outputBuffers.wrs->StoreKmer(kmerSeq, valuesToSave);

			if (outputBuffers.wrsSignificant)
			{
				if (pValue <= maxCorrectedPval)
				{
					outputBuffers.wrsSignificant->StoreKmer(kmerSeq, valuesToSave);
					outputBuffers.wrsSignificantCntMatrix->StoreKmer(kmerSeq, original_counts);
					outputBuffers.wrsSignificantFasta->StoreKmer(kmerSeq, binId, kmerIdInBin);
				}
			}
		}
		if (outputBuffers.dids)
		{
			const auto value = outEntry[valuesIdx++];
			outputBuffers.dids->StoreKmer(kmerSeq, value);
			keepNLargestCollection->addDids(kmerSeq, kmer, value, original_counts);
		}
		if (outputBuffers.anova)
		{
			const auto pValue = outEntry[valuesIdx++];
			const auto statistic = outEntry[valuesIdx++];

			const std::vector<Statistics_T> valuesToSave{ pValue, statistic };

			outputBuffers.anova->StoreKmer(kmerSeq, valuesToSave);

			if (outputBuffers.anovaSignificant)
			{
				if (pValue <= maxCorrectedPval)
				{
					outputBuffers.anovaSignificant->StoreKmer(kmerSeq, valuesToSave);
					outputBuffers.anovaSignificantCntMatrix->StoreKmer(kmerSeq, original_counts);
					outputBuffers.anovaSignificantFasta->StoreKmer(kmerSeq, binId, kmerIdInBin);
				}
			}
		}
	}
}


template<typename Statistics_T, typename VALUE_T>
void WritingGatherer<Statistics_T, VALUE_T>::initWriting(
	const std::unique_ptr<kmcdb::MetadataReader>& matrixMetadataReader,
	const std::vector<std::string>& cnt_matrix_output_header)
{
	const bool multiThreadedGeneration = params.mkmcParams.nThreads > 1;

	if (statisticsToGeneration.pearson)
	{
		writers.pearson = std::make_unique<TextFileWriter>(params.mkmcParams.outputFilePearson, multiThreadedGeneration);
		writers.pearson->StoreHeader({ "pearson" });
		Logger::Inst().Log("Info: generating Pearson correlation to " + params.mkmcParams.outputFilePearson + ".", 2);
		Logger::Inst().Log("Info: the file contains k-mers with the correlation values.", 2);
	}
	if (statisticsToGeneration.spearman)
	{
		writers.spearman = std::make_unique<TextFileWriter>(params.mkmcParams.outputFileSpearman, multiThreadedGeneration);
		writers.spearman->StoreHeader({ "spearman" });
		Logger::Inst().Log("Info: generating Spearman correlation to " + params.mkmcParams.outputFileSpearman + ".", 2);
		Logger::Inst().Log("Info: the file contains k-mers with the correlation values.", 2);
	}
	if (statisticsToGeneration.kendall)
	{
		writers.kendall = std::make_unique<TextFileWriter>(params.mkmcParams.outputFileKendall, multiThreadedGeneration);
		writers.kendall->StoreHeader({ "kendall" });
		Logger::Inst().Log("Info: generating Kendall Tau correlation to " + params.mkmcParams.outputFileKendall + ".", 2);
		Logger::Inst().Log("Info: the file contains k-mers with the correlation values.", 2);
	}
	if (statisticsToGeneration.entropy)
	{
		writers.entropy = std::make_unique<TextFileWriter>(params.mkmcParams.outputFileEntropy, multiThreadedGeneration);
		writers.entropy->StoreHeader({ "entropy" });
		Logger::Inst().Log("Info: generating entropy to " + params.mkmcParams.outputFileEntropy + ".", 2);
		Logger::Inst().Log("Info: the file contains k-mers with the entropy values.", 2);
	}
	if (statisticsToGeneration.tTest)
	{
		if (params.statisticsParams.correctPvalues)
		{
			writers.tTest = std::make_unique<TextFileWriter>(params.mkmcParams.outputFileTTestCor, multiThreadedGeneration);
			writers.tTest->StoreHeader({ "ttest_p_val_corrected", "ttest_df", "ttest_statistic" });
			Logger::Inst().Log("Info: generating corrected T-Test results to " + params.mkmcParams.outputFileTTestCor + ".", 2);
			Logger::Inst().Log("Info: the file contains k-mers with T-Test corrected p-values, df statistics, and T-Test values.", 2);

			writers.tTestSignificant = std::make_unique<TextFileWriter>(params.mkmcParams.outputFileTTestCorSignificant, multiThreadedGeneration);
			writers.tTestSignificant->StoreHeader({ "ttest_p_val_corrected", "ttest_df", "ttest_statistic" });
			Logger::Inst().Log("Info: generating corrected T-Test results to " + params.mkmcParams.outputFileTTestCorSignificant + ".", 2);
			Logger::Inst().Log("Info: the file contains statistically significant k-mers with T-Test corrected p-values, df statistics, and T-Test values.", 2);

			writers.tTestSignificantCntMatrix = std::make_unique<TextFileWriter>(params.mkmcParams.outputFileTTestCorSignificantCntMatrix, multiThreadedGeneration);
			writers.tTestSignificantCntMatrix->StoreHeader(cnt_matrix_output_header);
			Logger::Inst().Log("Info: generating matrix of statistically significant k-mers obtainted by corrected T-Test results to " + params.mkmcParams.outputFileTTestCorSignificantCntMatrix + ".", 2);
			Logger::Inst().Log("Info: the file contains counts matrix of statistically significant k-mers.", 2);

			writers.tTestSignificantFasta = std::make_unique<TextFileWriter>(params.mkmcParams.outputFileTTestCorSignificantFasta, multiThreadedGeneration);
			Logger::Inst().Log("Info: generating k-mers obtainted by corrected T-Test results to " + params.mkmcParams.outputFileTTestCorSignificantFasta + ".", 2);
			Logger::Inst().Log("Info: the file contains statistically significant k-mers in FASTA format.", 2);
		}
		else
		{
			writers.tTest = std::make_unique<TextFileWriter>(params.mkmcParams.outputFileTTest, multiThreadedGeneration);
			writers.tTest->StoreHeader({ "ttest_p_val", "ttest_df", "ttest_statistic" });
			Logger::Inst().Log("Info: generating T-Test results to " + params.mkmcParams.outputFileTTest + ".", 2);
			Logger::Inst().Log("Info: the file contains k-mers with T-Test p-values, df statistics, and T-Test values.", 2);
		}
	}
	if (statisticsToGeneration.snr)
	{
		writers.snr = std::make_unique<TextFileWriter>(params.mkmcParams.outputFileSNR, multiThreadedGeneration);
		writers.snr->StoreHeader({ "snr" });
		Logger::Inst().Log("Info: generating signal to noise ratio results to " + params.mkmcParams.outputFileSNR + ".", 2);
		Logger::Inst().Log("Info: the file contains k-mers with SNR results.", 2);
	}
	if (statisticsToGeneration.unnormalizedSnr)
	{
		writers.unnormalizedSnr = std::make_unique<TextFileWriter>(params.mkmcParams.outputFileUnnormalizedSNR, multiThreadedGeneration);
		writers.unnormalizedSnr->StoreHeader({ "snr_for_unnormalized" });
		Logger::Inst().Log("Info: generating unnormalized signal to noise ratio results to " + params.mkmcParams.outputFileUnnormalizedSNR + ".", 2);
		Logger::Inst().Log("Info: the file contains k-mers with SNR results.", 2);
	}
	if (statisticsToGeneration.wrs)
	{
		if (params.statisticsParams.correctPvalues)
		{
			writers.wrs = std::make_unique<TextFileWriter>(params.mkmcParams.outputFileWilcoxonRankSumCor, multiThreadedGeneration);
			writers.wrs->StoreHeader({ "wrs_p_val_corrected", "wrs_U1_statistic", "wrs_U2_statistic" });
			Logger::Inst().Log("Info: generating corrected Wilcoxon-rank sum (Mann-Whitney U test) results to " + params.mkmcParams.outputFileWilcoxonRankSumCor + ".", 2);
			Logger::Inst().Log("Info: the file contains k-mers with WRS corrected p-values, U1, and U2 statistics.", 2);

			writers.wrsSignificant = std::make_unique<TextFileWriter>(params.mkmcParams.outputFileWilcoxonRankSumCorSignificant, multiThreadedGeneration);
			writers.wrsSignificant->StoreHeader({ "wrs_p_val_corrected", "wrs_U1_statistic", "wrs_U2_statistic" });
			Logger::Inst().Log("Info: generating corrected Wilcoxon-rank sum (Mann-Whitney U test) results to " + params.mkmcParams.outputFileWilcoxonRankSumCorSignificant + ".", 2);
			Logger::Inst().Log("Info: the file contains statistically significant k-mers with WRS corrected p-values, U1, and U2 statistics.", 2);

			writers.wrsSignificantCntMatrix = std::make_unique<TextFileWriter>(params.mkmcParams.outputFileWilcoxonRankSumCorSignificantCntMatrix, multiThreadedGeneration);
			writers.wrsSignificantCntMatrix->StoreHeader(cnt_matrix_output_header);
			Logger::Inst().Log("Info: generating matrix of k-mers obtainted by corrected Wilcoxon-rank sum (Mann-Whitney U test) results to " + params.mkmcParams.outputFileWilcoxonRankSumCorSignificantCntMatrix + ".", 2);
			Logger::Inst().Log("Info: the file contains counts matrix of statistically significant k-mers.", 2);

			writers.wrsSignificantFasta = std::make_unique<TextFileWriter>(params.mkmcParams.outputFileWilcoxonRankSumCorSignificantFasta, multiThreadedGeneration);
			Logger::Inst().Log("Info: generating k-mers obtainted by corrected Wilcoxon-rank sum (Mann-Whitney U test) results to " + params.mkmcParams.outputFileWilcoxonRankSumCorSignificantFasta + ".", 2);
			Logger::Inst().Log("Info: the file contains statistically significant k-mers in FASTA format.", 2);
		}
		else
		{
			writers.wrs = std::make_unique<TextFileWriter>(params.mkmcParams.outputFileWilcoxonRankSum, multiThreadedGeneration);
			writers.wrs->StoreHeader({ "wrs_p_val", "wrs_U1_statistic", "wrs_U2_statistic" });
			Logger::Inst().Log("Info: generating Wilcoxon-rank sum (Mann-Whitney U test) results to " + params.mkmcParams.outputFileWilcoxonRankSum + ".", 2);
			Logger::Inst().Log("Info: the file contains k-mers with WRS p-values, U1, and U2 statistics.", 2);
		}
	}
	if (statisticsToGeneration.dids)
	{
		writers.dids = std::make_unique<TextFileWriter>(params.mkmcParams.outputFileDIDS, multiThreadedGeneration);
		writers.dids->StoreHeader({ "dids" });
		Logger::Inst().Log("Info: generating DIDS results to " + params.mkmcParams.outputFileDIDS + ".", 2);
		Logger::Inst().Log("Info: the file contains k-mers with DIDS results.", 2);
	}
	if (statisticsToGeneration.anova)
	{
		if (params.statisticsParams.correctPvalues)
		{
			writers.anova = std::make_unique<TextFileWriter>(params.mkmcParams.outputFileANOVACor, multiThreadedGeneration);
			writers.anova->StoreHeader({ "anova_p_val_corrected", "anova_statistic" });
			Logger::Inst().Log("Info: generating corrected ANOVA results to " + params.mkmcParams.outputFileANOVACor + ".", 2);
			Logger::Inst().Log("Info: the file contains k-mers with ANOVA p-values and statistic results.", 2);

			writers.anovaSignificant = std::make_unique<TextFileWriter>(params.mkmcParams.outputFileANOVACorSignificant, multiThreadedGeneration);
			writers.anovaSignificant->StoreHeader({ "anova_p_val_corrected", "anova_statistic" });
			Logger::Inst().Log("Info: generating corrected ANOVA results to " + params.mkmcParams.outputFileANOVACorSignificant + ".", 2);
			Logger::Inst().Log("Info: the file contains statistically significant k-mers with ANOVA corrected p-values and statistic results.", 2);

			writers.anovaSignificantCntMatrix = std::make_unique<TextFileWriter>(params.mkmcParams.outputFileANOVACorSignificantCntMatrix, multiThreadedGeneration);
			writers.anovaSignificantCntMatrix->StoreHeader(cnt_matrix_output_header);
			Logger::Inst().Log("Info: generating matrix of k-mers obtainted by corrected ANOVA results to " + params.mkmcParams.outputFileANOVACorSignificantCntMatrix + ".", 2);
			Logger::Inst().Log("Info: the file contains counts matrix of statistically significant k-mers.", 2);

			writers.anovaSignificantFasta = std::make_unique<TextFileWriter>(params.mkmcParams.outputFileANOVACorSignificantFasta, multiThreadedGeneration);
			Logger::Inst().Log("Info: generating k-mers obtainted by corrected ANOVA results to " + params.mkmcParams.outputFileANOVACorSignificantFasta + ".", 2);
			Logger::Inst().Log("Info: the file contains statistically significant k-mers in FASTA format.", 2);
		}
		else
		{
			writers.anova = std::make_unique<TextFileWriter>(params.mkmcParams.outputFileANOVA, multiThreadedGeneration);
			writers.anova->StoreHeader({ "anova_p_val", "anova_statistic" });
			Logger::Inst().Log("Info: generating ANOVA results to " + params.mkmcParams.outputFileANOVA + ".", 2);
			Logger::Inst().Log("Info: the file contains top k-mers with ANOVA p-values and statistic results.", 2);
		}
	}
}
