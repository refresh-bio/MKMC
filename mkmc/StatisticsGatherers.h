#pragma once

#include <memory>
#include <vector>
#include <string>
#include "parameters.h"
#include "DumpWriter.h"
#include "kmcdb/kmcdb.h"



struct StatisticsToGeneration
{
	bool normalize = false;

	bool pearson = false;
	bool spearman = false;
	bool kendall = false;

	bool entropy = false;

	bool differentialAnalysis = false; // logical sum of the following ones

	bool tTest = false;
	bool snr = false;
	bool wilcoxonRankSum = false;

	bool dids = false;
	bool anova = false;

	uint32_t nStatistics = 0;
	uint32_t nStatisticsWithPValues = 0;
	uint32_t nResults = 0; // sum of nStatistics and number of normalized samples
};



template<typename Statistics_T>
class WritingGatherer;



template<typename Statistics_T>
class WritingGathererBin
{
	friend WritingGatherer<Statistics_T>;
	WritingGatherer<Statistics_T>& mainWritingGatherer;

	const uint32_t kmerLength;
	const uint32_t numSamples;
	const double maxCorrectedPval;

	std::unique_ptr<OutputBuffer> pearsonOutputBuffer, spearmanOutputBuffer, kendallOutputBuffer;
	std::unique_ptr<OutputBuffer> entropyOutputBuffer;

	std::unique_ptr<OutputBuffer> tTestOutputBuffer;
	std::unique_ptr<OutputBuffer> tTestSignificantOutputBuffer;
	std::unique_ptr<OutputBuffer> tTestSignificantCntMatrixOutputBuffer;
	std::unique_ptr<OutputBuffer> tTestSignificantFastaOutputBuffer;

	std::unique_ptr<OutputBuffer> snrOutputBuffer;

	std::unique_ptr<OutputBuffer> wilcoxonRankSumOutputBuffer;
	std::unique_ptr<OutputBuffer> wilcoxonRankSumSignificantOutputBuffer;
	std::unique_ptr<OutputBuffer> wilcoxonRankSumSignificantCntMatrixOutputBuffer;
	std::unique_ptr<OutputBuffer> wilcoxonRankSumSignificantFastaOutputBuffer;

	std::unique_ptr<OutputBuffer> didsOutputBuffer;

	std::unique_ptr<OutputBuffer> anovaOutputBuffer;
	std::unique_ptr<OutputBuffer> anovaSignificantOutputBuffer;
	std::unique_ptr<OutputBuffer> anovaSignificantCntMatrixOutputBuffer;
	std::unique_ptr<OutputBuffer> anovaSignificantFastaOutputBuffer;

	kmcdb::BinWriterSortedPlain<Statistics_T>* outBin;

	const StatisticsToGeneration& statisticsToGeneration;

	size_t getMaxLineLength() const
	{
		//return kmerLength + 1 + refresh::numeric_conversion_max_length<Statistics_T>();
		return kmerLength + 1 + refresh::numeric_conversion_max_length<Statistics_T>() + 1; //@Maciej: I have added 1 since we also store EOL after value, right?
	}

	size_t getMaxLineLengthForFasta() const
	{
		//     >\n   k-mer       \n
		return 2   + kmerLength + 1;
	}
	size_t getMaxLineLengthForCntMatrix() const
	{
		constexpr uint32_t assumed_max_len_for_cnt = 20;
		//     k-mer      term(\t)            cnt                       term(\t or \n)
		return kmerLength + 1 + numSamples * (assumed_max_len_for_cnt + 1);
	}

	WritingGathererBin(WritingGatherer<Statistics_T>& mainWritingGatherer,
		const StatisticsToGeneration& statisticsToGeneration,
		const uint32_t kmerLength,
		uint32_t binId,
		uint32_t numSamples,
		double maxCorrectedPval);

public:
	template<unsigned SIZE, typename VALUE_T>
	void writeKmer(
		const std::vector<Statistics_T>& outEntry, // outEntry - normalized values (if any) followed by statistics
		const kmcdb::CKmer<SIZE>& kmer,
		const std::string kmerSeq,
		const std::vector<VALUE_T>& original_counts);
};



template<typename Statistics_T>
class WritingGatherer
{
	friend WritingGathererBin<Statistics_T>;
	Params& params;

	struct
	{
		std::unique_ptr<DumpWriter> pearson;
		std::unique_ptr<DumpWriter> spearman;
		std::unique_ptr<DumpWriter> kendall;

		std::unique_ptr<DumpWriter> entropy;

		std::unique_ptr<DumpWriter> tTest;
		std::unique_ptr<DumpWriter> tTestSignificant;
		std::unique_ptr<DumpWriter> tTestSignificantCntMatrix;
		std::unique_ptr<DumpWriter> tTestSignificantFasta;

		std::unique_ptr<DumpWriter> snr;
		std::unique_ptr<DumpWriter> wilcoxonRankSum;
		std::unique_ptr<DumpWriter> wilcoxonRankSumSignificant;
		std::unique_ptr<DumpWriter> wilcoxonRankSumSignificantCntMatrix;
		std::unique_ptr<DumpWriter> wilcoxonRankSumSignificantFasta;

		std::unique_ptr<DumpWriter> dids;
		std::unique_ptr<DumpWriter> anova;
		std::unique_ptr<DumpWriter> anovaSignificant;
		std::unique_ptr<DumpWriter> anovaSignificantCntMatrix;
		std::unique_ptr<DumpWriter> anovaSignificantFasta;
	} writers;

	std::unique_ptr<kmcdb::WriterSortedPlain<Statistics_T>> kmcdbWriter;

	const StatisticsToGeneration& statisticsToGeneration;

public:
	WritingGatherer(Params& params, const StatisticsToGeneration& statisticsToGeneration) :
		params(params),
		statisticsToGeneration(statisticsToGeneration)
	{}

	void initWriting(const std::unique_ptr<kmcdb::MetadataReader>& matrixMetadataReader, std::vector<std::string> sample_names); // pass sample_names by value

	std::unique_ptr<WritingGathererBin<Statistics_T>> getBin(uint32_t binId)
	{
		// explicit new calling is necessary due to constructor's privacy
		return std::unique_ptr<WritingGathererBin<Statistics_T>>(new WritingGathererBin<Statistics_T>(
			*this,
			statisticsToGeneration,
			params.stage1Params.GetKmerLen(),
			binId,
			static_cast<uint32_t>(params.mkmcParams.samples.size()),
			params.statisticsParams.maxCorrectedPval
		));
	}
};



template<typename Statistics_T>
WritingGathererBin<Statistics_T>::WritingGathererBin(WritingGatherer<Statistics_T>& mainWritingGatherer,
	const StatisticsToGeneration& statisticsToGeneration,
	const uint32_t kmerLength,
	uint32_t binId,
	uint32_t numSamples,
	double maxCorrectedPval) :
	mainWritingGatherer(mainWritingGatherer),
	kmerLength(kmerLength),
	numSamples(numSamples),
	maxCorrectedPval(maxCorrectedPval),
	outBin(mainWritingGatherer.kmcdbWriter->GetBin(binId)),
	statisticsToGeneration(statisticsToGeneration)
{
	if (statisticsToGeneration.pearson)
		pearsonOutputBuffer = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.pearson, getMaxLineLength());
	if (statisticsToGeneration.spearman)
		spearmanOutputBuffer = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.spearman, getMaxLineLength());
	if (statisticsToGeneration.kendall)
		kendallOutputBuffer = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.kendall, getMaxLineLength());

	if (statisticsToGeneration.entropy)
		entropyOutputBuffer = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.entropy, getMaxLineLength());

	if (statisticsToGeneration.tTest)
	{
		tTestOutputBuffer = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.tTest, getMaxLineLength());
		if (mainWritingGatherer.writers.tTestSignificant)
		{
			assert(mainWritingGatherer.writers.tTestSignificantCntMatrix);
			assert(mainWritingGatherer.writers.tTestSignificantFasta);

			tTestSignificantOutputBuffer = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.tTestSignificant, getMaxLineLength());
			tTestSignificantCntMatrixOutputBuffer = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.tTestSignificantCntMatrix, getMaxLineLengthForCntMatrix());
			tTestSignificantFastaOutputBuffer = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.tTestSignificantFasta, getMaxLineLengthForFasta());
		}
	}
	if (statisticsToGeneration.snr)
		snrOutputBuffer = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.snr, getMaxLineLength());
	if (statisticsToGeneration.wilcoxonRankSum)
	{
		wilcoxonRankSumOutputBuffer = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.wilcoxonRankSum, getMaxLineLength());
		if (mainWritingGatherer.writers.wilcoxonRankSumSignificant)
		{
			assert(mainWritingGatherer.writers.wilcoxonRankSumSignificantCntMatrix);
			assert(mainWritingGatherer.writers.wilcoxonRankSumSignificantFasta);

			wilcoxonRankSumSignificantOutputBuffer = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.wilcoxonRankSumSignificant, getMaxLineLength());
			wilcoxonRankSumSignificantCntMatrixOutputBuffer = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.wilcoxonRankSumSignificantCntMatrix, getMaxLineLengthForCntMatrix());
			wilcoxonRankSumSignificantFastaOutputBuffer = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.wilcoxonRankSumSignificantFasta, getMaxLineLengthForFasta());
		}
	}

	if (statisticsToGeneration.dids)
		didsOutputBuffer = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.dids, getMaxLineLength());
	if (statisticsToGeneration.anova)
	{
		anovaOutputBuffer = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.anova, getMaxLineLength());
		if (mainWritingGatherer.writers.anovaSignificant)
		{
			assert(mainWritingGatherer.writers.anovaSignificantCntMatrix);
			assert(mainWritingGatherer.writers.anovaSignificantFasta);

			anovaSignificantOutputBuffer = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.anovaSignificant, getMaxLineLength());
			anovaSignificantCntMatrixOutputBuffer = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.anovaSignificantCntMatrix, getMaxLineLengthForCntMatrix());
			anovaSignificantFastaOutputBuffer = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.anovaSignificantFasta, getMaxLineLengthForFasta());
		}
	}
}



template<typename Statistics_T>
template<unsigned SIZE, typename VALUE_T>
void WritingGathererBin<Statistics_T>::writeKmer(
	const std::vector<Statistics_T>& outEntry,
	const kmcdb::CKmer<SIZE>& kmer,
	const std::string kmerSeq,
	const std::vector<VALUE_T>& original_counts)
{
	size_t valuesIdx = statisticsToGeneration.nResults - statisticsToGeneration.nStatistics;
	if (mainWritingGatherer.statisticsToGeneration.pearson)
	{
		pearsonOutputBuffer->StoreKmer(kmerSeq, outEntry[valuesIdx++], StoreMethods::AsMatrixRow_single_val);
	}
	if (mainWritingGatherer.statisticsToGeneration.spearman)
	{
		spearmanOutputBuffer->StoreKmer(kmerSeq, outEntry[valuesIdx++], StoreMethods::AsMatrixRow_single_val);
	}
	if (mainWritingGatherer.statisticsToGeneration.kendall)
	{
		kendallOutputBuffer->StoreKmer(kmerSeq, outEntry[valuesIdx++], StoreMethods::AsMatrixRow_single_val);
	}

	if (mainWritingGatherer.statisticsToGeneration.entropy)
	{
		entropyOutputBuffer->StoreKmer(kmerSeq, outEntry[valuesIdx++], StoreMethods::AsMatrixRow_single_val);
	}
	if (mainWritingGatherer.statisticsToGeneration.differentialAnalysis)
	{
		if (mainWritingGatherer.statisticsToGeneration.tTest)
		{
			auto value = outEntry[valuesIdx++];
			tTestOutputBuffer->StoreKmer(kmerSeq, value, StoreMethods::AsMatrixRow_single_val);

			if (tTestSignificantOutputBuffer)
			{
				if (value <= maxCorrectedPval)
				{
					tTestSignificantOutputBuffer->StoreKmer(kmerSeq, value, StoreMethods::AsMatrixRow_single_val);
					tTestSignificantCntMatrixOutputBuffer->StoreKmer(kmerSeq, original_counts, StoreMethods::AsMatrixRow);
					tTestSignificantFastaOutputBuffer->StoreKmer(kmerSeq, original_counts, StoreMethods::AsFastaRecord);
				}
			}
		}
		if (mainWritingGatherer.statisticsToGeneration.snr)
		{
			snrOutputBuffer->StoreKmer(kmerSeq, outEntry[valuesIdx++], StoreMethods::AsMatrixRow_single_val);
		}
		if (mainWritingGatherer.statisticsToGeneration.wilcoxonRankSum)
		{
			auto value = outEntry[valuesIdx++];
			wilcoxonRankSumOutputBuffer->StoreKmer(kmerSeq, value, StoreMethods::AsMatrixRow_single_val);

			if (wilcoxonRankSumSignificantOutputBuffer)
			{
				if (value <= maxCorrectedPval)
				{
					wilcoxonRankSumOutputBuffer->StoreKmer(kmerSeq, value, StoreMethods::AsMatrixRow_single_val);
					wilcoxonRankSumSignificantCntMatrixOutputBuffer->StoreKmer(kmerSeq, original_counts, StoreMethods::AsMatrixRow);
					wilcoxonRankSumSignificantFastaOutputBuffer->StoreKmer(kmerSeq, original_counts, StoreMethods::AsFastaRecord);
				}
			}
		}
		if (mainWritingGatherer.statisticsToGeneration.dids)
		{
			didsOutputBuffer->StoreKmer(kmerSeq, outEntry[valuesIdx++], StoreMethods::AsMatrixRow_single_val);
		}
		if (mainWritingGatherer.statisticsToGeneration.anova)
		{
			auto value = outEntry[valuesIdx++];
			anovaOutputBuffer->StoreKmer(kmerSeq, value, StoreMethods::AsMatrixRow_single_val);

			if (anovaSignificantOutputBuffer)
			{
				if (value <= maxCorrectedPval)
				{
					anovaSignificantOutputBuffer->StoreKmer(kmerSeq, value, StoreMethods::AsMatrixRow_single_val);
					anovaSignificantCntMatrixOutputBuffer->StoreKmer(kmerSeq, original_counts, StoreMethods::AsMatrixRow);
					anovaSignificantFastaOutputBuffer->StoreKmer(kmerSeq, original_counts, StoreMethods::AsFastaRecord);
				}
			}
		}
	}

	outBin->AddKmer(kmer, outEntry.data());
}



template<typename Statistics_T>
void WritingGatherer<Statistics_T>::initWriting(const std::unique_ptr<kmcdb::MetadataReader>& matrixMetadataReader, std::vector<std::string> sample_names)
{
	kmcdb::Config config;
	config.num_bins = matrixMetadataReader->GetConfig().num_bins;
	config.signature_len = matrixMetadataReader->GetConfig().signature_len;
	config.signature_selection_scheme = matrixMetadataReader->GetConfig().signature_selection_scheme;
	config.signature_to_bin_mapping = matrixMetadataReader->GetConfig().signature_to_bin_mapping;
	config.kmer_len = matrixMetadataReader->GetConfig().kmer_len;
	config.num_samples = sample_names.size();
	config.num_bytes_single_value = { sizeof(Statistics_T) };

	config.num_samples += params.statisticsParams.correlationMethods.size(); //I will add this correlations as a new columns
	config.num_samples += params.statisticsParams.generateEntropy ? 1 : 0;
	config.num_samples += params.statisticsParams.classificationMethods.size();
	//this make sense because those all of the same type (currently double)

	const bool multiThreadedGeneration = params.mkmcParams.nThreads > 1;

	std::vector<std::string> cnt_matrix_output_header = sample_names;

	if (statisticsToGeneration.pearson)
	{
		writers.pearson = std::make_unique<DumpWriter>(params.mkmcParams.outputFilePearson, multiThreadedGeneration);
		writers.pearson->StoreHeader({ "pearson_cor" });
		sample_names.emplace_back("pearson_cor");
	}
	if (statisticsToGeneration.spearman)
	{
		writers.spearman = std::make_unique<DumpWriter>(params.mkmcParams.outputFileSpearman, multiThreadedGeneration);
		writers.spearman->StoreHeader({ "spearman_cor" });
		sample_names.emplace_back("spearman_cor");
	}
	if (statisticsToGeneration.kendall)
	{
		writers.kendall = std::make_unique<DumpWriter>(params.mkmcParams.outputFileKendall, multiThreadedGeneration);
		writers.kendall->StoreHeader({ "kendall_cor" });
		sample_names.emplace_back("kendall_cor");
	}
	if (statisticsToGeneration.entropy)
	{
		writers.entropy = std::make_unique<DumpWriter>(params.mkmcParams.outputFileEntropy, multiThreadedGeneration);
		writers.entropy->StoreHeader({ "entropy" });
		sample_names.emplace_back("entropy");
	}
	if (statisticsToGeneration.tTest)
	{
		if (params.statisticsParams.correctPvalues)
		{
			writers.tTest = std::make_unique<DumpWriter>(params.mkmcParams.outputFileTTestCor, multiThreadedGeneration);
			writers.tTest->StoreHeader({ "ttest_analysis_p_val_cor" });
			sample_names.emplace_back("ttest_analysis_cor");

			writers.tTestSignificant = std::make_unique<DumpWriter>(params.mkmcParams.outputFileTTestCorSignificant, multiThreadedGeneration);
			writers.tTestSignificant->StoreHeader({ "ttest_analysis_p_val_cor" });

			writers.tTestSignificantCntMatrix = std::make_unique<DumpWriter>(params.mkmcParams.outputFileTTestCorSignificantCntMatrix, multiThreadedGeneration);
			writers.tTestSignificantCntMatrix->StoreHeader(cnt_matrix_output_header);

			writers.tTestSignificantFasta = std::make_unique<DumpWriter>(params.mkmcParams.outputFileTTestCorSignificantFasta, multiThreadedGeneration);

		}
		else
		{
			writers.tTest = std::make_unique<DumpWriter>(params.mkmcParams.outputFileTTest, multiThreadedGeneration);
			writers.tTest->StoreHeader({ "ttest_analysis_p_val" });
			sample_names.emplace_back("ttest_analysis");
		}
	}
	if (statisticsToGeneration.snr)
	{
		writers.snr = std::make_unique<DumpWriter>(params.mkmcParams.outputFileSNR, multiThreadedGeneration);
		writers.snr->StoreHeader({ "snr_analysis" });
		sample_names.emplace_back("snr_analysis");
	}
	if (statisticsToGeneration.wilcoxonRankSum)
	{
		if (params.statisticsParams.correctPvalues)
		{
			writers.wilcoxonRankSum = std::make_unique<DumpWriter>(params.mkmcParams.outputFileWilcoxonRankSumCor, multiThreadedGeneration);
			writers.wilcoxonRankSum->StoreHeader({ "wrs_analysis_p_val_cor" });
			sample_names.emplace_back("wrs_analysis_cor");

			writers.wilcoxonRankSumSignificant = std::make_unique<DumpWriter>(params.mkmcParams.outputFileWilcoxonRankSumCorSignificant, multiThreadedGeneration);
			writers.wilcoxonRankSumSignificant->StoreHeader({ "wrs_analysis_p_val_cor" });

			writers.wilcoxonRankSumSignificantCntMatrix = std::make_unique<DumpWriter>(params.mkmcParams.outputFileWilcoxonRankSumCorSignificantCntMatrix, multiThreadedGeneration);
			writers.wilcoxonRankSumSignificantCntMatrix->StoreHeader(cnt_matrix_output_header);

			writers.wilcoxonRankSumSignificantFasta = std::make_unique<DumpWriter>(params.mkmcParams.outputFileWilcoxonRankSumCorSignificantFasta, multiThreadedGeneration);
		}
		else
		{
			writers.wilcoxonRankSum = std::make_unique<DumpWriter>(params.mkmcParams.outputFileWilcoxonRankSum, multiThreadedGeneration);
			writers.wilcoxonRankSum->StoreHeader({ "wrs_analysis_p_val" });
			sample_names.emplace_back("wrs_analysis");
		}
	}
	if (statisticsToGeneration.dids)
	{
		writers.dids = std::make_unique<DumpWriter>(params.mkmcParams.outputFileDIDS, multiThreadedGeneration);
		writers.dids->StoreHeader({ "dids_analysis" });
		sample_names.emplace_back("dids_analysis");
	}
	if (statisticsToGeneration.anova)
	{
		if (params.statisticsParams.correctPvalues)
		{
			writers.anova = std::make_unique<DumpWriter>(params.mkmcParams.outputFileANOVACor, multiThreadedGeneration);
			writers.anova->StoreHeader({ "anova_analysis_p_val_cor" });
			sample_names.emplace_back("anova_analysis_cor");

			writers.anovaSignificant = std::make_unique<DumpWriter>(params.mkmcParams.outputFileANOVACorSignificant, multiThreadedGeneration);
			writers.anovaSignificant->StoreHeader({ "anova_analysis_p_val_cor" });

			writers.anovaSignificantCntMatrix = std::make_unique<DumpWriter>(params.mkmcParams.outputFileANOVACorSignificantCntMatrix, multiThreadedGeneration);
			writers.anovaSignificantCntMatrix->StoreHeader(cnt_matrix_output_header);

			writers.anovaSignificantFasta = std::make_unique<DumpWriter>(params.mkmcParams.outputFileANOVACorSignificantFasta, multiThreadedGeneration);
		}
		else
		{
			writers.anova = std::make_unique<DumpWriter>(params.mkmcParams.outputFileANOVA, multiThreadedGeneration);
			writers.anova->StoreHeader({ "anova_analysis_p_val" });
			sample_names.emplace_back("anova_analysis");
		}
	}

	kmcdb::ConfigSortedPlain representation_config{};

	kmcdbWriter = std::make_unique<kmcdb::WriterSortedPlain<Statistics_T>>(
		config,
		representation_config,
		params.mkmcParams.outputStatsBinFile,
		params.mkmcParams.outputBinFile,
		sample_names);
}
