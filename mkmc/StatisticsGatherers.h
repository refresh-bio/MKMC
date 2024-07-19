#pragma once

#include <memory>
#include <vector>
#include <string>
#include "parameters.h"
#include "DumpWriter.h"
#include "kmcdb/kmcdb.h"



struct StatisticsToGeneration
{
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
};



template<typename Statistics_T>
class WritingGatherer;



template<typename Statistics_T>
class WritingGathererBin
{
	friend WritingGatherer<Statistics_T>;
	WritingGatherer<Statistics_T>& mainWritingGatherer;

	const uint32_t kmerLength;

	std::unique_ptr<OutputBuffer> pearsonOutputBuffer, spearmanOutputBuffer, kendallOutputBuffer;
	std::unique_ptr<OutputBuffer> entropyOutputBuffer;
	std::unique_ptr<OutputBuffer> tTestOutputBuffer, snrOutputBuffer, wilcoxonRankSumOutputBuffer;
	std::unique_ptr<OutputBuffer> didsOutputBuffer, anovaOutputBuffer;

	kmcdb::BinWriterSortedPlain<Statistics_T>* outBin;

	const StatisticsToGeneration& statisticsToGeneration;

	WritingGathererBin(WritingGatherer<Statistics_T>& mainWritingGatherer, const StatisticsToGeneration& statisticsToGeneration, const uint32_t kmerLength, uint32_t binId);

public:
	template<unsigned SIZE>
	void writeKmer(std::vector<Statistics_T>& outNormMatrixEntry, const kmcdb::CKmer<SIZE>& kmer, const std::string kmerSeq, const std::vector<Statistics_T>& outStatsEntry); // outNormMatrixEntry may be modified
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
		std::unique_ptr<DumpWriter> snr;
		std::unique_ptr<DumpWriter> wilcoxonRankSum;

		std::unique_ptr<DumpWriter> dids;
		std::unique_ptr<DumpWriter> anova;
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
		return std::unique_ptr<WritingGathererBin<Statistics_T>>(new WritingGathererBin<Statistics_T>(*this, statisticsToGeneration, params.stage1Params.GetKmerLen(), binId));
	}
};



template<typename Statistics_T>
WritingGathererBin<Statistics_T>::WritingGathererBin(WritingGatherer<Statistics_T>& mainWritingGatherer, const StatisticsToGeneration& statisticsToGeneration, const uint32_t kmerLength, uint32_t binId) :
	mainWritingGatherer(mainWritingGatherer),
	kmerLength(kmerLength),
	outBin(mainWritingGatherer.kmcdbWriter->GetBin(binId)),
	statisticsToGeneration(statisticsToGeneration)
{
	if (statisticsToGeneration.pearson)
		pearsonOutputBuffer = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.pearson, kmerLength);
	if (statisticsToGeneration.spearman)
		spearmanOutputBuffer = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.spearman, kmerLength);
	if (statisticsToGeneration.kendall)
		kendallOutputBuffer = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.kendall, kmerLength);

	if (statisticsToGeneration.entropy)
		entropyOutputBuffer = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.entropy, kmerLength);

	if (statisticsToGeneration.tTest)
		tTestOutputBuffer = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.tTest, kmerLength);
	if (statisticsToGeneration.snr)
		snrOutputBuffer = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.snr, kmerLength);
	if (statisticsToGeneration.wilcoxonRankSum)
		wilcoxonRankSumOutputBuffer = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.wilcoxonRankSum, kmerLength);

	if (statisticsToGeneration.dids)
		didsOutputBuffer = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.dids, kmerLength);
	if (statisticsToGeneration.anova)
		anovaOutputBuffer = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.anova, kmerLength);
}



template<typename Statistics_T>
template<unsigned SIZE>
void WritingGathererBin<Statistics_T>::writeKmer(std::vector<Statistics_T>& outNormMatrixEntry, const kmcdb::CKmer<SIZE>& kmer, const std::string kmerSeq, const std::vector<Statistics_T>& outStatsEntry)
{
	auto storeMethod = []<typename VALUE_T>(const std::string & kmerSeq, const VALUE_T cnt, char* out) -> size_t
	{
		std::memcpy(out, kmerSeq.data(), kmerSeq.length());
		out += kmerSeq.length();
		*out = '\t';
		++out;

		size_t res = kmerSeq.length() + 1;

		auto store_single_value = [&](const VALUE_T& val, char term)
		{
			size_t r{};
			if constexpr (std::is_integral_v<VALUE_T>)
				r = refresh::int_to_pchar(val, out, term);
			else if constexpr (std::is_floating_point_v<VALUE_T>)
			{
				if (std::isnan(val))
				{
					out[0] = 'n';
					out[1] = 'a';
					out[2] = 'n';
					out[3] = term;
					r = 4;
				}
				else if (std::isinf(val))
				{
					if (val < static_cast<VALUE_T>(0)) {
						out[0] = '-';
						out[1] = 'i';
						out[2] = 'n';
						out[3] = 'f';
						out[4] = term;
						r = 5;
					}
					else
					{
						out[0] = 'i';
						out[1] = 'n';
						out[2] = 'f';
						out[3] = term;
						r = 4;
					}
				}
				else
					r = refresh::real_to_pchar(val, out, 6, term);

			}
			else
			{
				static_assert(!sizeof(VALUE_T), "Unsupported type");
			}
			out += r;
			res += r;
		};

		store_single_value(cnt, '\n');

		return res;
	};

	size_t valuesIdx = 0;
	if (mainWritingGatherer.statisticsToGeneration.pearson)
	{
		outNormMatrixEntry.push_back(outStatsEntry[valuesIdx]);
		pearsonOutputBuffer->StoreKmer(kmerSeq, outStatsEntry[valuesIdx], storeMethod);
		++valuesIdx;
	}
	if (mainWritingGatherer.statisticsToGeneration.spearman)
	{
		outNormMatrixEntry.push_back(outStatsEntry[valuesIdx]);
		spearmanOutputBuffer->StoreKmer(kmerSeq, outStatsEntry[valuesIdx], storeMethod);
		++valuesIdx;
	}
	if (mainWritingGatherer.statisticsToGeneration.kendall)
	{
		outNormMatrixEntry.push_back(outStatsEntry[valuesIdx]);
		kendallOutputBuffer->StoreKmer(kmerSeq, outStatsEntry[valuesIdx], storeMethod);
		++valuesIdx;
	}

	if (mainWritingGatherer.statisticsToGeneration.entropy)
	{
		outNormMatrixEntry.push_back(outStatsEntry[valuesIdx]);
		entropyOutputBuffer->StoreKmer(kmerSeq, outStatsEntry[valuesIdx], storeMethod);
		++valuesIdx;
	}
	if (mainWritingGatherer.statisticsToGeneration.differentialAnalysis)
	{
		if (mainWritingGatherer.statisticsToGeneration.tTest)
		{
			outNormMatrixEntry.push_back(outStatsEntry[valuesIdx]);
			tTestOutputBuffer->StoreKmer(kmerSeq, outStatsEntry[valuesIdx], storeMethod);
			++valuesIdx;
		}
		if (mainWritingGatherer.statisticsToGeneration.snr)
		{
			outNormMatrixEntry.push_back(outStatsEntry[valuesIdx]);
			snrOutputBuffer->StoreKmer(kmerSeq, outStatsEntry[valuesIdx], storeMethod);
			++valuesIdx;
		}
		if (mainWritingGatherer.statisticsToGeneration.wilcoxonRankSum)
		{
			outNormMatrixEntry.push_back(outStatsEntry[valuesIdx]);
			wilcoxonRankSumOutputBuffer->StoreKmer(kmerSeq, outStatsEntry[valuesIdx], storeMethod);
			++valuesIdx;
		}
		if (mainWritingGatherer.statisticsToGeneration.dids)
		{
			outNormMatrixEntry.push_back(outStatsEntry[valuesIdx]);
			didsOutputBuffer->StoreKmer(kmerSeq, outStatsEntry[valuesIdx], storeMethod);
			++valuesIdx;
		}
		if (mainWritingGatherer.statisticsToGeneration.anova)
		{
			outNormMatrixEntry.push_back(outStatsEntry[valuesIdx]);
			anovaOutputBuffer->StoreKmer(kmerSeq, outStatsEntry[valuesIdx], storeMethod);
			++valuesIdx;
		}
	}

	outBin->AddKmer(kmer, outNormMatrixEntry.data());
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
	config.num_samples = matrixMetadataReader->GetConfig().num_samples;
	config.num_bytes_single_value = { sizeof(Statistics_T) };

	config.num_samples += params.statisticsParams.correlationMethods.size(); //I will add this correlations as a new columns
	config.num_samples += params.statisticsParams.generateEntropy ? 1 : 0;
	config.num_samples += params.statisticsParams.classificationMethods.size();
	//this make sense because those all of the same type (currently double)

	const bool multiThreadedGeneration = params.mkmcParams.nThreads > 1;

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
		writers.tTest = std::make_unique<DumpWriter>(params.mkmcParams.outputFileTTest, multiThreadedGeneration);
		writers.tTest->StoreHeader({ "ttest_analysis_p_val" });
		sample_names.emplace_back("ttest_analysis");
	}
	if (statisticsToGeneration.snr)
	{
		writers.snr = std::make_unique<DumpWriter>(params.mkmcParams.outputFileSNR, multiThreadedGeneration);
		writers.snr->StoreHeader({ "snr_analysis" });
		sample_names.emplace_back("snr_analysis");
	}
	if (statisticsToGeneration.wilcoxonRankSum)
	{
		writers.wilcoxonRankSum = std::make_unique<DumpWriter>(params.mkmcParams.outputFileWilcoxonRankSum, multiThreadedGeneration);
		writers.wilcoxonRankSum->StoreHeader({ "wrs_analysis_p_val" });
		sample_names.emplace_back("wrs_analysis");
	}
	if (statisticsToGeneration.dids)
	{
		writers.dids = std::make_unique<DumpWriter>(params.mkmcParams.outputFileDIDS, multiThreadedGeneration);
		writers.dids->StoreHeader({ "dids_analysis" });
		sample_names.emplace_back("dids_analysis");
	}
	if (statisticsToGeneration.anova)
	{
		writers.anova = std::make_unique<DumpWriter>(params.mkmcParams.outputFileANOVA, multiThreadedGeneration);
		writers.anova->StoreHeader({ "anova_analysis_p_val" });
		sample_names.emplace_back("anova_analysis");
	}

	kmcdb::ConfigSortedPlain representation_config{};

	kmcdbWriter = std::make_unique<kmcdb::WriterSortedPlain<Statistics_T>>(
		config,
		representation_config,
		params.mkmcParams.outputStatsBinFile,
		params.mkmcParams.outputBinFile,
		sample_names);
}
