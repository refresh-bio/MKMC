#pragma once

#include <memory>
#include <vector>
#include <string>

#include "Merger.h"
#include "parameters.h"
#include "DumpWriter.h"
#include "kmcdb/kmcdb.h"
#include "KeepNLargests.h"
#include "lib/refresh/statistics/lib/statistics_umap.h"


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
	bool unnormalizedSnr = false;
	bool wrs = false;

	bool dids = false;
	bool anova = false;

	uint32_t nStatistics = 0;
	uint32_t nStatisticsWithPValues = 0;
	uint32_t nAdditionalValuesOfCorrectedStats = 0; // T-Test = 2, WRS = 2, ANOVA = 1
};



template<typename Statistics_T>
class WritingGatherer;




//mkokot_TODO: move to some other file?

//mkokot_TODO: its not perfect because I am copying this vector with values
//This could be avoided if I copy it only if my element is large enough to
//be keep in the heap
//but it would require some extension in KeepNLargest probably
template<unsigned SIZE>
struct KeepTopElem
{
	std::string kmerSeq;
	kmcdb::CKmer<SIZE> kmer;
	double key;
	std::vector<uint64_t> counts;

	KeepTopElem(std::string kmerSeq,
	kmcdb::CKmer<SIZE> kmer,
	double key,
	std::vector<uint64_t> counts):
	kmerSeq(std::move(kmerSeq)),
	kmer(kmer),
	key(key),
	counts(std::move(counts))
	{

	}

	struct ABSGreater
	{
		bool operator()(const KeepTopElem<SIZE>& lhs, const KeepTopElem<SIZE>& rhs)
		{
			return std::make_pair(std::abs(lhs.key), rhs.kmer) > std::make_pair(std::abs(rhs.key), lhs.kmer);
		}
	};

	struct Greater
	{
		bool operator()(const KeepTopElem<SIZE>& lhs, const KeepTopElem<SIZE>& rhs)
		{
			return std::make_pair(lhs.key, rhs.kmer) > std::make_pair(rhs.key, lhs.kmer);
		}
	};
};




//mkokot_TODO: move it somewhere else?
template<unsigned SIZE>
struct KeepNLargestCollection
{
	using ABSGreater = typename KeepTopElem<SIZE>::ABSGreater;
	using Greater = typename KeepTopElem<SIZE>::Greater;

	using KeepTopNLargestABS_T = KeepNLargests<KeepTopElem<SIZE>, ABSGreater>;
	using KeepTopNLargestPlain_T = KeepNLargests<KeepTopElem<SIZE>, Greater>;

	std::unique_ptr<KeepTopNLargestABS_T> pearson;
	std::unique_ptr<KeepTopNLargestABS_T> spearman;
	std::unique_ptr<KeepTopNLargestABS_T> kendall;

	std::unique_ptr<KeepTopNLargestPlain_T> entropy;

	std::unique_ptr<KeepTopNLargestABS_T> snr;
	std::unique_ptr<KeepTopNLargestABS_T> unnormalizedSnr;
	std::unique_ptr<KeepTopNLargestPlain_T> dids;
};

template<unsigned SIZE>
struct KeepNLargestCollectionGlobal
{
private:
	std::mutex mtx;
	KeepNLargestCollection<SIZE> global;

	template<typename PRED>
	static void add_for(
		std::unique_ptr<KeepNLargests<KeepTopElem<SIZE>, PRED>>& src,
		std::unique_ptr<KeepNLargests<KeepTopElem<SIZE>, PRED>>& dest)
	{
		//if source was not collected do nothing
		if (!src)
			return;

		if (!dest)
			dest = std::make_unique<KeepNLargests<KeepTopElem<SIZE>, PRED>>(src->GetN());

		assert(dest->GetN() == src->GetN()); //just to be sure that all source have the same N

		std::vector<KeepTopElem<SIZE>> data;
		src->Steal(data);
		for (auto& elem : data)
			dest->Add(std::move(elem));
	}

	template<typename PRED>
	static void flush_for(
		std::unique_ptr<KeepNLargests<KeepTopElem<SIZE>, PRED>>& to_flush,
		const std::string& fname_top, const std::vector<std::string>& header_top, size_t max_line_len_top,
		const std::string& fname_top_matrix, const std::vector<std::string>& header_top_matrix, size_t max_line_len_top_matrix,
		const std::string& fname_top_fasta, size_t max_line_len_top_fasta)
	{
		if (!to_flush)
			return;

		DumpWriter writer_top(fname_top, false);
		writer_top.StoreHeader(header_top);
		OutputBuffer buff_top(writer_top, max_line_len_top);

		DumpWriter writer_top_matrix(fname_top_matrix, false);
		writer_top_matrix.StoreHeader(header_top_matrix);
		OutputBuffer buff_top_matrix(writer_top_matrix, max_line_len_top_matrix);

		DumpWriter writer_fasta(fname_top_fasta, false);
		OutputBuffer buff_top_fasta(writer_fasta, max_line_len_top_fasta);

		std::vector<KeepTopElem<SIZE>> data;
		to_flush->StealSorted(data, PRED{}); //could actually be Steal (no sorted), but lets keep it deterministic
		for (auto& elem : data)
		{
			buff_top.StoreKmer(elem.kmerSeq, elem.key, StoreMethods::AsMatrixRow_single_val);
			buff_top_matrix.StoreKmer(elem.kmerSeq, elem.counts, StoreMethods::AsMatrixRow);
			buff_top_fasta.StoreKmer(elem.kmerSeq, elem.counts, StoreMethods::AsFastaRecord);
		}
	}
public:
	void Add(KeepNLargestCollection<SIZE>& collection)
	{
		std::lock_guard lck(mtx);
		add_for(collection.pearson, global.pearson);
		add_for(collection.spearman, global.spearman);
		add_for(collection.kendall, global.kendall);

		add_for(collection.entropy, global.entropy);

		add_for(collection.snr, global.snr);
		add_for(collection.unnormalizedSnr, global.unnormalizedSnr);
		add_for(collection.dids, global.dids);
	}

	void Flush(const Params& params, const std::vector<std::string>& cnt_matrix_output_header)
	{
		//mkokot_TODO: ugly code repetition, we have the same (almost) in other class as private methods, to be refactored...
		auto getMaxLineLength = [](uint32_t kmerLength) -> size_t
		{
			return kmerLength + 1 + refresh::numeric_conversion_max_length<double>() + 1;
		};

		auto getMaxLineLengthForFasta = [](uint32_t kmerLength) -> size_t
		{
			//     >\n   k-mer       \n
			return 2 + kmerLength + 1;
		};
		auto getMaxLineLengthForCntMatrix = [](uint32_t kmerLength, uint32_t numSamples) -> size_t
		{
			constexpr uint32_t assumed_max_len_for_cnt = 20;
			//     k-mer      term(\t)            cnt                       term(\t or \n)
			return kmerLength + 1 + numSamples * (assumed_max_len_for_cnt + 1);
		};

		auto max_line_len_top = getMaxLineLength(params.stage1Params.GetKmerLen());
		auto max_line_len_top_matrix = getMaxLineLengthForCntMatrix(params.stage1Params.GetKmerLen(), static_cast<uint32_t>(params.mkmcParams.samples.size()));
		auto max_line_len_top_fasta = getMaxLineLengthForFasta(params.stage1Params.GetKmerLen());

		flush_for(global.pearson,
			params.mkmcParams.outputFilePearsonTop, { "pearson_cor" }, max_line_len_top,
			params.mkmcParams.outputFilePearsonTopCntMatrix, cnt_matrix_output_header, max_line_len_top_matrix,
			params.mkmcParams.outputFilePearsonTopFasta, max_line_len_top_fasta);

		flush_for(global.spearman,
			params.mkmcParams.outputFileSpearmanTop, { "spearman_cor" }, max_line_len_top,
			params.mkmcParams.outputFileSpearmanTopCntMatrix, cnt_matrix_output_header, max_line_len_top_matrix,
			params.mkmcParams.outputFileSpearmanTopFasta, max_line_len_top_fasta);

		flush_for(global.kendall,
			params.mkmcParams.outputFileKendallTop, { "kendall_cor" }, max_line_len_top,
			params.mkmcParams.outputFileKendallTopCntMatrix, cnt_matrix_output_header, max_line_len_top_matrix,
			params.mkmcParams.outputFileKendallTopFasta, max_line_len_top_fasta);

		flush_for(global.entropy,
			params.mkmcParams.outputFileEntropyTop, { "entropy" }, max_line_len_top,
			params.mkmcParams.outputFileEntropyTopCntMatrix, cnt_matrix_output_header, max_line_len_top_matrix,
			params.mkmcParams.outputFileEntropyTopFasta, max_line_len_top_fasta);

		flush_for(global.snr,
			params.mkmcParams.outputFileSNRTop, { "snr_analysis" }, max_line_len_top,
			params.mkmcParams.outputFileSNRTopCntMatrix, cnt_matrix_output_header, max_line_len_top_matrix,
			params.mkmcParams.outputFileSNRTopFasta, max_line_len_top_fasta);

		flush_for(global.unnormalizedSnr,
			params.mkmcParams.outputFileUnnormalizedSNRTop, { "snr_analysis_for_unnormalized" }, max_line_len_top,
			params.mkmcParams.outputFileUnnormalizedSNRTopCntMatrix, cnt_matrix_output_header, max_line_len_top_matrix,
			params.mkmcParams.outputFileUnnormalizedSNRTopFasta, max_line_len_top_fasta);

		flush_for(global.dids,
			params.mkmcParams.outputFileDIDSTop, { "dids_analysis" }, max_line_len_top,
			params.mkmcParams.outputFileDIDSTopCntMatrix, cnt_matrix_output_header, max_line_len_top_matrix,
			params.mkmcParams.outputFileDIDSTopFasta, max_line_len_top_fasta);
	}
};


template<typename Statistics_T>
class WritingGathererBin
{
	friend WritingGatherer<Statistics_T>;
	WritingGatherer<Statistics_T>& mainWritingGatherer;

	const uint32_t kmerLength;
	const uint32_t numSamples;
	const double maxCorrectedPval;

	struct {
		std::unique_ptr<OutputBuffer> pearson;
		std::unique_ptr<OutputBuffer> spearman;
		std::unique_ptr<OutputBuffer> kendall;

		std::unique_ptr<OutputBuffer> entropy;

		std::unique_ptr<OutputBuffer> tTest;
		std::unique_ptr<OutputBuffer> tTestSignificant;
		std::unique_ptr<OutputBuffer> tTestSignificantCntMatrix;
		std::unique_ptr<OutputBuffer> tTestSignificantFasta;

		std::unique_ptr<OutputBuffer> snr;
		std::unique_ptr<OutputBuffer> unnormalizedSnr;

		std::unique_ptr<OutputBuffer> wrs;
		std::unique_ptr<OutputBuffer> wrsSignificant;
		std::unique_ptr<OutputBuffer> wrsSignificantCntMatrix;
		std::unique_ptr<OutputBuffer> wrsSignificantFasta;

		std::unique_ptr<OutputBuffer> dids;

		std::unique_ptr<OutputBuffer> anova;
		std::unique_ptr<OutputBuffer> anovaSignificant;
		std::unique_ptr<OutputBuffer> anovaSignificantCntMatrix;
		std::unique_ptr<OutputBuffer> anovaSignificantFasta;
	} outputBuffers;

	size_t getMaxLineLength() const
	{
		return kmerLength + 1 + refresh::numeric_conversion_max_length<Statistics_T>() + 1; //+ 1 since we also store EOL after value
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
		const uint32_t kmerLength,
		uint32_t binId,
		uint32_t numSamples,
		double maxCorrectedPval,
		bool gatherCorrected,
		bool gatherNotCorrected);

public:
	template<unsigned SIZE, typename VALUE_T>
	void writeKmer(
		const std::vector<Statistics_T>& outEntry, // outEntry - normalized values (if any) followed by statistics
		const kmcdb::CKmer<SIZE>& kmer,
		const std::string& kmerSeq,
		const std::vector<VALUE_T>& original_counts,
		KeepNLargestCollection<SIZE>* keepNLargestCollection = nullptr);

};



template<typename Statistics_T>
class WritingGatherer
{
	friend WritingGathererBin<Statistics_T>;
	const Params& params;

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
		std::unique_ptr<DumpWriter> unnormalizedSnr;

		std::unique_ptr<DumpWriter> wrs;
		std::unique_ptr<DumpWriter> wrsSignificant;
		std::unique_ptr<DumpWriter> wrsSignificantCntMatrix;
		std::unique_ptr<DumpWriter> wrsSignificantFasta;

		std::unique_ptr<DumpWriter> dids;

		std::unique_ptr<DumpWriter> anova;
		std::unique_ptr<DumpWriter> anovaSignificant;
		std::unique_ptr<DumpWriter> anovaSignificantCntMatrix;
		std::unique_ptr<DumpWriter> anovaSignificantFasta;
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

	std::unique_ptr<WritingGathererBin<Statistics_T>> getBin(
		uint32_t binId,
		bool gatherCorrected = true,
		bool gatherNotCorrected = true)
	{
		// explicit new calling is necessary due to constructor's privacy
		return std::unique_ptr<WritingGathererBin<Statistics_T>>(new WritingGathererBin<Statistics_T>(
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



template<typename Statistics_T>
WritingGathererBin<Statistics_T>::WritingGathererBin(
	WritingGatherer<Statistics_T>& mainWritingGatherer,
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
		outputBuffers.pearson = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.pearson, getMaxLineLength());
	if (gatherNotCorrected && mainWritingGatherer.writers.spearman)
		outputBuffers.spearman = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.spearman, getMaxLineLength());
	if (gatherNotCorrected && mainWritingGatherer.writers.kendall)
		outputBuffers.kendall = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.kendall, getMaxLineLength());

	if (gatherNotCorrected && mainWritingGatherer.writers.entropy)
		outputBuffers.entropy = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.entropy, getMaxLineLength());

	if (gatherCorrected && mainWritingGatherer.writers.tTest)
	{
		outputBuffers.tTest = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.tTest, getMaxLineLength());
		if (mainWritingGatherer.writers.tTestSignificant)
		{
			assert(mainWritingGatherer.writers.tTestSignificantCntMatrix);
			assert(mainWritingGatherer.writers.tTestSignificantFasta);

			outputBuffers.tTestSignificant = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.tTestSignificant, getMaxLineLength());
			outputBuffers.tTestSignificantCntMatrix = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.tTestSignificantCntMatrix, getMaxLineLengthForCntMatrix());
			outputBuffers.tTestSignificantFasta = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.tTestSignificantFasta, getMaxLineLengthForFasta());
		}
	}
	if (gatherNotCorrected && mainWritingGatherer.writers.snr)
		outputBuffers.snr = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.snr, getMaxLineLength());
	if (gatherNotCorrected && mainWritingGatherer.writers.unnormalizedSnr)
		outputBuffers.unnormalizedSnr = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.unnormalizedSnr, getMaxLineLength());
	if (gatherCorrected && mainWritingGatherer.writers.wrs)
	{
		outputBuffers.wrs = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.wrs, getMaxLineLength());
		if (mainWritingGatherer.writers.wrsSignificant)
		{
			assert(mainWritingGatherer.writers.wrsSignificantCntMatrix);
			assert(mainWritingGatherer.writers.wrsSignificantFasta);

			outputBuffers.wrsSignificant = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.wrsSignificant, getMaxLineLength());
			outputBuffers.wrsSignificantCntMatrix = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.wrsSignificantCntMatrix, getMaxLineLengthForCntMatrix());
			outputBuffers.wrsSignificantFasta = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.wrsSignificantFasta, getMaxLineLengthForFasta());
		}
	}

	if (gatherNotCorrected && mainWritingGatherer.writers.dids)
		outputBuffers.dids = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.dids, getMaxLineLength());
	if (gatherCorrected && mainWritingGatherer.writers.anova)
	{
		outputBuffers.anova = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.anova, getMaxLineLength());
		if (mainWritingGatherer.writers.anovaSignificant)
		{
			assert(mainWritingGatherer.writers.anovaSignificantCntMatrix);
			assert(mainWritingGatherer.writers.anovaSignificantFasta);

			outputBuffers.anovaSignificant = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.anovaSignificant, getMaxLineLength());
			outputBuffers.anovaSignificantCntMatrix = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.anovaSignificantCntMatrix, getMaxLineLengthForCntMatrix());
			outputBuffers.anovaSignificantFasta = std::make_unique<OutputBuffer>(*mainWritingGatherer.writers.anovaSignificantFasta, getMaxLineLengthForFasta());
		}
	}
}



template<typename Statistics_T>
template<unsigned SIZE, typename VALUE_T>
void WritingGathererBin<Statistics_T>::writeKmer(
	const std::vector<Statistics_T>& outEntry,
	const kmcdb::CKmer<SIZE>& kmer,
	const std::string& kmerSeq,
	const std::vector<VALUE_T>& original_counts,
	KeepNLargestCollection<SIZE>* keepNLargestCollection/* = nullptr*/)
{
	const bool safeNTop = keepNLargestCollection && mainWritingGatherer.params.statisticsParams.nTop;
	size_t valuesIdx = 0;
	if (outputBuffers.pearson)
	{
		const auto value = outEntry[valuesIdx++];
		outputBuffers.pearson->StoreKmer(kmerSeq, value, StoreMethods::AsMatrixRow_single_val);

		if (safeNTop)
			keepNLargestCollection->pearson->Add(KeepTopElem<SIZE>{kmerSeq, kmer, value, original_counts});
	}
	if (outputBuffers.spearman)
	{
		const auto value = outEntry[valuesIdx++];
		outputBuffers.spearman->StoreKmer(kmerSeq, value, StoreMethods::AsMatrixRow_single_val);

		if (safeNTop)
			keepNLargestCollection->spearman->Add(KeepTopElem<SIZE>{kmerSeq, kmer, value, original_counts});
	}
	if (outputBuffers.kendall)
	{
		const auto value = outEntry[valuesIdx++];
		outputBuffers.kendall->StoreKmer(kmerSeq, value, StoreMethods::AsMatrixRow_single_val);

		if (safeNTop)
			keepNLargestCollection->kendall->Add(KeepTopElem<SIZE>{kmerSeq, kmer, value, original_counts});
	}

	if (outputBuffers.entropy)
	{
		const auto value = outEntry[valuesIdx++];
		outputBuffers.entropy->StoreKmer(kmerSeq, value, StoreMethods::AsMatrixRow_single_val);

		if (safeNTop)
			keepNLargestCollection->entropy->Add(KeepTopElem<SIZE>{kmerSeq, kmer, value, original_counts});
	}
	if (mainWritingGatherer.statisticsToGeneration.differentialAnalysis)
	{
		if (outputBuffers.tTest)
		{
			const auto pValue = outEntry[valuesIdx++];
			const auto df = outEntry[valuesIdx++];
			const auto statistic = outEntry[valuesIdx++];

			const std::vector<double> valuesToSave{ pValue, df, statistic };

			outputBuffers.tTest->StoreKmer(kmerSeq, valuesToSave, StoreMethods::AsMatrixRow);

			if (outputBuffers.tTestSignificant)
			{
				if (pValue <= maxCorrectedPval)
				{
					outputBuffers.tTestSignificant->StoreKmer(kmerSeq, valuesToSave, StoreMethods::AsMatrixRow);
					outputBuffers.tTestSignificantCntMatrix->StoreKmer(kmerSeq, original_counts, StoreMethods::AsMatrixRow);
					outputBuffers.tTestSignificantFasta->StoreKmer(kmerSeq, original_counts, StoreMethods::AsFastaRecord);
				}
			}
		}
		if (outputBuffers.snr)
		{
			const auto value = outEntry[valuesIdx++];
			outputBuffers.snr->StoreKmer(kmerSeq, value, StoreMethods::AsMatrixRow_single_val);

			if (safeNTop)
				keepNLargestCollection->snr->Add(KeepTopElem<SIZE>{kmerSeq, kmer, value, original_counts});
		}
		if (outputBuffers.unnormalizedSnr)
		{
			const auto value = outEntry[valuesIdx++];
			outputBuffers.unnormalizedSnr->StoreKmer(kmerSeq, value, StoreMethods::AsMatrixRow_single_val);

			if (safeNTop)
				keepNLargestCollection->unnormalizedSnr->Add(KeepTopElem<SIZE>{kmerSeq, kmer, value, original_counts});
		}
		if (outputBuffers.wrs)
		{
			const auto pValue = outEntry[valuesIdx++];
			const auto statisticU1 = outEntry[valuesIdx++];
			const auto statisticU2 = outEntry[valuesIdx++];

			const std::vector<double> valuesToSave{ pValue, statisticU1, statisticU2 };

			outputBuffers.wrs->StoreKmer(kmerSeq, valuesToSave, StoreMethods::AsMatrixRow);

			if (outputBuffers.wrsSignificant)
			{
				if (pValue <= maxCorrectedPval)
				{
					outputBuffers.wrsSignificant->StoreKmer(kmerSeq, valuesToSave, StoreMethods::AsMatrixRow);
					outputBuffers.wrsSignificantCntMatrix->StoreKmer(kmerSeq, original_counts, StoreMethods::AsMatrixRow);
					outputBuffers.wrsSignificantFasta->StoreKmer(kmerSeq, original_counts, StoreMethods::AsFastaRecord);
				}
			}
		}
		if (outputBuffers.dids)
		{
			const auto value = outEntry[valuesIdx++];
			outputBuffers.dids->StoreKmer(kmerSeq, value, StoreMethods::AsMatrixRow_single_val);

			if (safeNTop)
				keepNLargestCollection->dids->Add(KeepTopElem<SIZE>{kmerSeq, kmer, value, original_counts});
		}
		if (outputBuffers.anova)
		{
			const auto pValue = outEntry[valuesIdx++];
			const auto statistic = outEntry[valuesIdx++];

			const std::vector<double> valuesToSave{ pValue, statistic };

			outputBuffers.anova->StoreKmer(kmerSeq, valuesToSave, StoreMethods::AsMatrixRow);

			if (outputBuffers.anovaSignificant)
			{
				if (pValue <= maxCorrectedPval)
				{
					outputBuffers.anovaSignificant->StoreKmer(kmerSeq, valuesToSave, StoreMethods::AsMatrixRow);
					outputBuffers.anovaSignificantCntMatrix->StoreKmer(kmerSeq, original_counts, StoreMethods::AsMatrixRow);
					outputBuffers.anovaSignificantFasta->StoreKmer(kmerSeq, original_counts, StoreMethods::AsFastaRecord);
				}
			}
		}
	}
}


template<typename Statistics_T>
void WritingGatherer<Statistics_T>::initWriting(
	const std::unique_ptr<kmcdb::MetadataReader>& matrixMetadataReader,
	const std::vector<std::string>& cnt_matrix_output_header)
{
	const bool multiThreadedGeneration = params.mkmcParams.nThreads > 1;

	if (statisticsToGeneration.pearson)
	{
		writers.pearson = std::make_unique<DumpWriter>(params.mkmcParams.outputFilePearson, multiThreadedGeneration);
		writers.pearson->StoreHeader({ "pearson_cor" });
	}
	if (statisticsToGeneration.spearman)
	{
		writers.spearman = std::make_unique<DumpWriter>(params.mkmcParams.outputFileSpearman, multiThreadedGeneration);
		writers.spearman->StoreHeader({ "spearman_cor" });
	}
	if (statisticsToGeneration.kendall)
	{
		writers.kendall = std::make_unique<DumpWriter>(params.mkmcParams.outputFileKendall, multiThreadedGeneration);
		writers.kendall->StoreHeader({ "kendall_cor" });
	}
	if (statisticsToGeneration.entropy)
	{
		writers.entropy = std::make_unique<DumpWriter>(params.mkmcParams.outputFileEntropy, multiThreadedGeneration);
		writers.entropy->StoreHeader({ "entropy" });
	}
	if (statisticsToGeneration.tTest)
	{
		if (params.statisticsParams.correctPvalues)
		{
			writers.tTest = std::make_unique<DumpWriter>(params.mkmcParams.outputFileTTestCor, multiThreadedGeneration);
			writers.tTest->StoreHeader({ "ttest_analysis_p_val_cor", "ttest_analysis_df", "ttest_analysis_statistic" });

			writers.tTestSignificant = std::make_unique<DumpWriter>(params.mkmcParams.outputFileTTestCorSignificant, multiThreadedGeneration);
			writers.tTestSignificant->StoreHeader({ "ttest_analysis_p_val_cor", "ttest_analysis_df", "ttest_analysis_statistic" });

			writers.tTestSignificantCntMatrix = std::make_unique<DumpWriter>(params.mkmcParams.outputFileTTestCorSignificantCntMatrix, multiThreadedGeneration);
			writers.tTestSignificantCntMatrix->StoreHeader(cnt_matrix_output_header);

			writers.tTestSignificantFasta = std::make_unique<DumpWriter>(params.mkmcParams.outputFileTTestCorSignificantFasta, multiThreadedGeneration);
		}
		else
		{
			writers.tTest = std::make_unique<DumpWriter>(params.mkmcParams.outputFileTTest, multiThreadedGeneration);
			writers.tTest->StoreHeader({ "ttest_analysis_p_val", "ttest_analysis_df", "ttest_analysis_statistic" });
		}
	}
	if (statisticsToGeneration.snr)
	{
		writers.snr = std::make_unique<DumpWriter>(params.mkmcParams.outputFileSNR, multiThreadedGeneration);
		writers.snr->StoreHeader({ "snr_analysis" });
	}
	if (statisticsToGeneration.unnormalizedSnr)
	{
		writers.unnormalizedSnr = std::make_unique<DumpWriter>(params.mkmcParams.outputFileUnnormalizedSNR, multiThreadedGeneration);
		writers.unnormalizedSnr->StoreHeader({ "snr_analysis_for_unnormalized" });
	}
	if (statisticsToGeneration.wrs)
	{
		if (params.statisticsParams.correctPvalues)
		{
			writers.wrs = std::make_unique<DumpWriter>(params.mkmcParams.outputFileWilcoxonRankSumCor, multiThreadedGeneration);
			writers.wrs->StoreHeader({ "wrs_analysis_p_val_cor", "wrs_analysis_U1_statistic", "wrs_analysis_U2_statistic" });

			writers.wrsSignificant = std::make_unique<DumpWriter>(params.mkmcParams.outputFileWilcoxonRankSumCorSignificant, multiThreadedGeneration);
			writers.wrsSignificant->StoreHeader({ "wrs_analysis_p_val_cor", "wrs_analysis_U1_statistic", "wrs_analysis_U2_statistic" });

			writers.wrsSignificantCntMatrix = std::make_unique<DumpWriter>(params.mkmcParams.outputFileWilcoxonRankSumCorSignificantCntMatrix, multiThreadedGeneration);
			writers.wrsSignificantCntMatrix->StoreHeader(cnt_matrix_output_header);

			writers.wrsSignificantFasta = std::make_unique<DumpWriter>(params.mkmcParams.outputFileWilcoxonRankSumCorSignificantFasta, multiThreadedGeneration);
		}
		else
		{
			writers.wrs = std::make_unique<DumpWriter>(params.mkmcParams.outputFileWilcoxonRankSum, multiThreadedGeneration);
			writers.wrs->StoreHeader({ "wrs_analysis_p_val", "wrs_analysis_U1_statistic", "wrs_analysis_U2_statistic" });
		}
	}
	if (statisticsToGeneration.dids)
	{
		writers.dids = std::make_unique<DumpWriter>(params.mkmcParams.outputFileDIDS, multiThreadedGeneration);
		writers.dids->StoreHeader({ "dids_analysis" });
	}
	if (statisticsToGeneration.anova)
	{
		if (params.statisticsParams.correctPvalues)
		{
			writers.anova = std::make_unique<DumpWriter>(params.mkmcParams.outputFileANOVACor, multiThreadedGeneration);
			writers.anova->StoreHeader({ "anova_analysis_p_val_cor", "anova_statistic" });

			writers.anovaSignificant = std::make_unique<DumpWriter>(params.mkmcParams.outputFileANOVACorSignificant, multiThreadedGeneration);
			writers.anovaSignificant->StoreHeader({ "anova_analysis_p_val_cor", "anova_statistic" });

			writers.anovaSignificantCntMatrix = std::make_unique<DumpWriter>(params.mkmcParams.outputFileANOVACorSignificantCntMatrix, multiThreadedGeneration);
			writers.anovaSignificantCntMatrix->StoreHeader(cnt_matrix_output_header);

			writers.anovaSignificantFasta = std::make_unique<DumpWriter>(params.mkmcParams.outputFileANOVACorSignificantFasta, multiThreadedGeneration);
		}
		else
		{
			writers.anova = std::make_unique<DumpWriter>(params.mkmcParams.outputFileANOVA, multiThreadedGeneration);
			writers.anova->StoreHeader({ "anova_analysis_p_val", "anova_statistic" });
		}
	}

	kmcdb::ConfigSortedPlain representation_config{};
}
