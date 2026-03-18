#pragma once

#include <vector>
#include <algorithm>

template<typename T, typename PRED = std::greater<T>>
class KeepNLargests {
	std::vector<T> heap;
	size_t n;
	PRED pred;
public:
	KeepNLargests(size_t n, PRED pred = PRED{}) : n(n), pred(pred) {
		assert(n);
	}
	void Add(T&& elem) {
		if (heap.size() == n) {
			if (pred(heap[0], elem))
				return;
			heap.push_back(std::move(elem));
			std::push_heap(heap.begin(), heap.end(), pred);
			std::pop_heap(heap.begin(), heap.end(), pred);
			heap.pop_back();
		}
		else {
			heap.push_back(std::move(elem));
			if (heap.size() == n) {
				std::make_heap(heap.begin(), heap.end(), pred);
			}
		}
	}
	size_t GetN() const
	{
		return n;
	}
	const std::vector<T>& Get() const {
		return heap;
	}

	void Steal(std::vector<T>& res) const {
		res = std::move(heap);
	}

	template<typename Pred = std::less<T>>
	std::vector<T> GetSorted(Pred pred = std::less<T>{}) {
		std::vector<T> res = heap;
		std::sort(res.begin(), res.end(), pred);
		return res;
	}

	template<typename Pred = std::less<T>>
	void StealSorted(std::vector<T>& res, Pred pred = std::less<T>{}) {
		std::sort(heap.begin(), heap.end(), pred);
		res = std::move(heap);
	}
};



//mkokot_TODO: its not perfect because I am copying this vector with values
//This could be avoided if I copy it only if my element is large enough to
//be keep in the heap
//but it would require some extension in KeepNLargest probably
template<unsigned SIZE, typename Statistics_T, typename VALUE_T>
struct KeepTopElem
{
	// kmerSeq and counts are stored in the heap, however they are not necessary to determine proper order, but are useful for storing to files in the future
	std::string kmerSeq;
	kmcdb::CKmer<SIZE> kmer;
	Statistics_T key;
	std::vector<VALUE_T> counts;

	KeepTopElem(std::string kmerSeq,
		kmcdb::CKmer<SIZE> kmer,
		Statistics_T key,
		std::vector<VALUE_T> counts) :
		kmerSeq(std::move(kmerSeq)),
		kmer(kmer),
		key(key),
		counts(std::move(counts))
	{}

	struct ABSGreater
	{
		bool operator()(const KeepTopElem& lhs, const KeepTopElem& rhs)
		{
			if (std::isnan(lhs.key))
			{
				if (std::isnan(rhs.key))
					return lhs.kmer.seeded_murmur64() < rhs.kmer.seeded_murmur64();
				return false;
			}
			if (std::isnan(rhs.key))
				return true;
			return std::make_pair(std::abs(lhs.key), rhs.kmer.seeded_murmur64()) > std::make_pair(std::abs(rhs.key), lhs.kmer.seeded_murmur64());
		}
	};

	struct Greater
	{
		bool operator()(const KeepTopElem& lhs, const KeepTopElem& rhs)
		{
			if (std::isnan(lhs.key))
			{
				if (std::isnan(rhs.key))
					return lhs.kmer.seeded_murmur64() < rhs.kmer.seeded_murmur64();
				return false;
			}
			if (std::isnan(rhs.key))
				return true;
			return std::make_pair(lhs.key, rhs.kmer.seeded_murmur64()) > std::make_pair(rhs.key, lhs.kmer.seeded_murmur64());
		}
	};
};



template<unsigned SIZE, typename Statistics_T, typename VALUE_T>
class KeepNLargestCollectionBase
{
protected:
	using Elem = KeepTopElem<SIZE, Statistics_T, VALUE_T>;

	using ABSGreater = typename Elem::ABSGreater;
	using Greater = typename Elem::Greater;

	using KeepTopNLargestABS_T = KeepNLargests<Elem, ABSGreater>;
	using KeepTopNLargestPlain_T = KeepNLargests<Elem, Greater>;

	template<typename PRED>
	static void add_for(
		KeepNLargests<Elem, PRED>& src,
		KeepNLargests<Elem, PRED>& dest);

	template<typename PRED>
	static void flush_for(
		KeepNLargests<Elem, PRED>& to_flush,
		uint32_t first_col_len, size_t num_columns,
		const std::string& fname_top, const std::vector<std::string>& header_top,
		const std::string& fname_top_matrix, const std::vector<std::string>& header_top_matrix,
		const std::string& fname_top_fasta);
};



template<unsigned SIZE, typename Statistics_T, typename VALUE_T>
class KeepNLargestCollection : KeepNLargestCollectionBase<SIZE, Statistics_T, VALUE_T>
{
	template<unsigned SIZE_, typename Statistics_T_, typename VALUE_T_, typename KeepNLargestCollection_T_>
	friend class KeepNLargestCollectionGlobal;

	using KeepNLargestCollectionBase<SIZE, Statistics_T, VALUE_T>::add_for;
	using KeepNLargestCollectionBase<SIZE, Statistics_T, VALUE_T>::flush_for;

	using typename KeepNLargestCollectionBase<SIZE, Statistics_T, VALUE_T>::Elem;
	using typename KeepNLargestCollectionBase<SIZE, Statistics_T, VALUE_T>::KeepTopNLargestABS_T;
	using typename KeepNLargestCollectionBase<SIZE, Statistics_T, VALUE_T>::KeepTopNLargestPlain_T;

	std::unique_ptr<KeepTopNLargestABS_T> pearson;
	std::unique_ptr<KeepTopNLargestABS_T> spearman;
	std::unique_ptr<KeepTopNLargestABS_T> kendall;

	std::unique_ptr<KeepTopNLargestPlain_T> entropy;

	std::unique_ptr<KeepTopNLargestPlain_T> snr;
	std::unique_ptr<KeepTopNLargestPlain_T> unnormalizedSnr;
	std::unique_ptr<KeepTopNLargestPlain_T> dids;

	KeepNLargestCollection()
	{}

	template<typename PRED>
	void add_impl(std::unique_ptr<KeepNLargests<Elem, PRED>>& src, std::unique_ptr<KeepNLargests<Elem, PRED>>& dest) const
	{
		//if source was not collected do nothing
		if (!src)
			return;

		if (!dest)
			dest = std::make_unique<KeepNLargests<Elem, PRED>>(src->GetN());

		add_for(*src, *dest);
	}

	void add(KeepNLargestCollection<SIZE, Statistics_T, VALUE_T>& collection)
	{
		add_impl(collection.pearson, pearson);
		add_impl(collection.spearman, spearman);
		add_impl(collection.kendall, kendall);

		add_impl(collection.entropy, entropy);

		add_impl(collection.snr, snr);
		add_impl(collection.unnormalizedSnr, unnormalizedSnr);
		add_impl(collection.dids, dids);
	}

	void flush(const Params& params, const std::vector<std::string>& cnt_matrix_output_header)
	{
		const size_t nSamples = params.mkmcParams.samples.size();
		const size_t kmerLen = params.stage1Params.GetKmerLen();

		bool tooLittleKmers = false;
		size_t nTopKmers = 0, nKmersExists = 0;

		auto logText = [](const std::string& alg, auto& top, auto& topCntMatrix, auto& topFasta)
		{
			Logger::Inst().Log("Info: generating " + alg + " top results to " + top + ".", 2);
			Logger::Inst().Log("Info: the file contains top k-mers with correlation values.", 2);
			Logger::Inst().Log("Info: generating " + alg + " top results to " + topCntMatrix + ".", 2);
			Logger::Inst().Log("Info: the file contains counts matrix of top k-mers.", 2);
			Logger::Inst().Log("Info: generating " + alg + " top results to " + topFasta + ".", 2);
			Logger::Inst().Log("Info: the file contains top k-mers in FASTA format.", 2);
		};

		auto flush = [nSamples, kmerLen, cnt_matrix_output_header, &tooLittleKmers, &nTopKmers, &nKmersExists, logText](const std::string& algHeader, const std::string& algorithmName, const auto& dataToFlush, auto& top, auto& topCntMatrix, auto& topFasta)
		{
			if (!dataToFlush)
				return;

			if (dataToFlush->GetN() > dataToFlush->Get().size())
			{
				assert((!tooLittleKmers && nKmersExists == 0) || nTopKmers == dataToFlush->GetN());
				assert((!tooLittleKmers && nTopKmers == 0) || nKmersExists == dataToFlush->Get().size());
				tooLittleKmers = true;
				
				nTopKmers = dataToFlush->GetN();
				nKmersExists = dataToFlush->Get().size();
			}

			flush_for(*dataToFlush, kmerLen, nSamples, top, { algHeader }, topCntMatrix, cnt_matrix_output_header, topFasta);
			logText(algorithmName, top, topCntMatrix, topFasta);
		};

		flush("pearson", "Pearson", pearson, params.mkmcParams.outputFilePearsonTop, params.mkmcParams.outputFilePearsonTopCntMatrix, params.mkmcParams.outputFilePearsonTopFasta);
		flush("spearman", "Spearman", spearman, params.mkmcParams.outputFileSpearmanTop, params.mkmcParams.outputFileSpearmanTopCntMatrix, params.mkmcParams.outputFileSpearmanTopFasta);
		flush("kendall", "Kendall Tau", kendall, params.mkmcParams.outputFileKendallTop, params.mkmcParams.outputFileKendallTopCntMatrix, params.mkmcParams.outputFileKendallTopFasta);
		flush("entropy", "entropy", entropy, params.mkmcParams.outputFileEntropyTop, params.mkmcParams.outputFileEntropyTopCntMatrix, params.mkmcParams.outputFileEntropyTopFasta);
		flush("snr", "signal to noise ratio", snr, params.mkmcParams.outputFileSNRTop, params.mkmcParams.outputFileSNRTopCntMatrix, params.mkmcParams.outputFileSNRTopFasta);
		flush("snr_for_unnormalized", "unnormalized signal to noise ratio", unnormalizedSnr, params.mkmcParams.outputFileUnnormalizedSNRTop, params.mkmcParams.outputFileUnnormalizedSNRTopCntMatrix, params.mkmcParams.outputFileUnnormalizedSNRTopFasta);
		flush("dids", "DIDS", dids, params.mkmcParams.outputFileDIDSTop, params.mkmcParams.outputFileDIDSTopCntMatrix, params.mkmcParams.outputFileDIDSTopFasta);

		if (tooLittleKmers)
			Logger::Inst().Log("Warning: requested to store " + std::to_string(nTopKmers) + " k-mers with top p-values, but only " + std::to_string(nKmersExists) + " k-mers exist.");
	}

public:
	KeepNLargestCollection(size_t nTop, bool pearson, bool spearman, bool kendall, bool entropy, bool snr, bool unnormalizedSnr, bool dids)
	{
		if (nTop == 0)
			return;

		if (pearson)
			this->pearson = std::make_unique<KeepTopNLargestABS_T>(nTop);

		if (spearman)
			this->spearman = std::make_unique<KeepTopNLargestABS_T>(nTop);

		if (kendall)
			this->kendall = std::make_unique<KeepTopNLargestABS_T>(nTop);

		if (entropy)
			this->entropy = std::make_unique<KeepTopNLargestPlain_T>(nTop);

		if (snr)
			this->snr = std::make_unique<KeepTopNLargestPlain_T>(nTop);

		if (unnormalizedSnr)
			this->unnormalizedSnr = std::make_unique<KeepTopNLargestPlain_T>(nTop);

		if (dids)
			this->dids = std::make_unique<KeepTopNLargestPlain_T>(nTop);
	}

	void addPearson(
		const std::string& kmerSeq,
		const kmcdb::CKmer<SIZE>& kmer,
		Statistics_T key,
		const std::vector<VALUE_T>& counts)
	{
		if (pearson)
			pearson->Add(Elem{ kmerSeq, kmer, key, counts });
	}

	void addSpearman(
		const std::string& kmerSeq,
		const kmcdb::CKmer<SIZE>& kmer,
		Statistics_T key,
		const std::vector<VALUE_T>& counts)
	{
		if (spearman)
			spearman->Add(Elem{ kmerSeq, kmer, key, counts });
	}

	void addKendall(
		const std::string& kmerSeq,
		const kmcdb::CKmer<SIZE>& kmer,
		Statistics_T key,
		const std::vector<VALUE_T>& counts)
	{
		if (kendall)
			kendall->Add(Elem{ kmerSeq, kmer, key, counts });
	}

	void addEntropy(
		const std::string& kmerSeq,
		const kmcdb::CKmer<SIZE>& kmer,
		Statistics_T key,
		const std::vector<VALUE_T>& counts)
	{
		if (entropy)
			entropy->Add(Elem{ kmerSeq, kmer, key, counts });
	}

	void addSnr(
		const std::string& kmerSeq,
		const kmcdb::CKmer<SIZE>& kmer,
		Statistics_T key,
		const std::vector<VALUE_T>& counts)
	{
		if (snr)
			snr->Add(Elem{ kmerSeq, kmer, key, counts });
	}

	void addUnnormalizedSnr(
		const std::string& kmerSeq,
		const kmcdb::CKmer<SIZE>& kmer,
		Statistics_T key,
		const std::vector<VALUE_T>& counts)
	{
		if (unnormalizedSnr)
			unnormalizedSnr->Add(Elem{ kmerSeq, kmer, key, counts });
	}

	void addDids(
		const std::string& kmerSeq,
		const kmcdb::CKmer<SIZE>& kmer,
		Statistics_T key,
		const std::vector<VALUE_T>& counts)
	{
		if (dids)
			dids->Add(Elem{ kmerSeq, kmer, key, counts });
	}
};



template<unsigned SIZE, typename Statistics_T, typename VALUE_T>
class KeepNLargestCollectionCV : KeepNLargestCollectionBase<SIZE, Statistics_T, VALUE_T>
{
	template<unsigned SIZE_, typename Statistics_T_, typename VALUE_T_, typename KeepNLargestCollection_T_>
	friend class KeepNLargestCollectionGlobal;

	using KeepNLargestCollectionBase<SIZE, Statistics_T, VALUE_T>::add_for;
	using KeepNLargestCollectionBase<SIZE, Statistics_T, VALUE_T>::flush_for;

	using typename KeepNLargestCollectionBase<SIZE, Statistics_T, VALUE_T>::Elem;
	using typename KeepNLargestCollectionBase<SIZE, Statistics_T, VALUE_T>::KeepTopNLargestABS_T;

	const std::vector<size_t>* const samplesToBeTestOrder;

	std::vector<std::unique_ptr<KeepTopNLargestABS_T>> pearson;
	std::vector<std::unique_ptr<KeepTopNLargestABS_T>> spearman;
	std::vector<std::unique_ptr<KeepTopNLargestABS_T>> kendall;

	KeepNLargestCollectionCV() : samplesToBeTestOrder(nullptr){};

	template<typename PRED>
	void add_impl(std::vector<std::unique_ptr<KeepNLargests<Elem, PRED>>>& src, std::vector<std::unique_ptr<KeepNLargests<Elem, PRED>>>& dest) const
	{
		if (src.empty())
			return;
		
		if (dest.empty())
		{
			dest.resize(src.size());
			for (auto& it : dest)
				it = std::make_unique<KeepNLargests<Elem, PRED>>(src.front()->GetN());
		}
		else
			assert(src.size() == dest.size());

		for (size_t i = 0; i < src.size(); ++i)
			add_for(*src[i], *dest[i]);
	}

	void add(KeepNLargestCollectionCV<SIZE, Statistics_T, VALUE_T>& collection)
	{
		add_impl(collection.pearson, pearson);
		add_impl(collection.spearman, spearman);
		add_impl(collection.kendall, kendall);
	}

	void flush(const Params& params, const std::vector<std::string>& whole_cnt_matrix_output_header)
	{
		if (!params.statisticsParams.cvParams.cv)
			return;

		const size_t nTestSamples = params.statisticsParams.cvParams.nTestSamples;
		const size_t nSamples = params.mkmcParams.samples.size();
		const size_t nTrainSamples = nSamples - nTestSamples;

		assert(nSamples % nTestSamples == 0);
		const size_t nCols = nSamples / nTestSamples;

		size_t nFolds = pearson.size();
		if (spearman.size() > 0)
		{
			assert(nFolds == 0 || nFolds == spearman.size());
			nFolds = spearman.size();
		}
		if (kendall.size() > 0)
		{
			assert(nFolds == 0 || nFolds == kendall.size());
			nFolds = kendall.size();
		}

		bool tooLittleKmers = false;
		size_t nTopKmers = 0, nKmersExists = 0;

		using enum StatisticsParams::CorrelationMethod;
		const std::vector<size_t>& samplesToBeTestOrder = params.statisticsParams.cvParams.samplesToBeTestOrder;

		std::vector<std::string> matrixHeader(nSamples); // method of matrixHeader content generation is similar as in CVGenerator

		size_t iTestSamples = 0;
		for (size_t iFold = 0; iFold < nFolds; ++iFold)
		{
			std::string trainSamplesText, testSamplesText;
			size_t iCurrentFoldTestSamples = 0;
			for (size_t iSamples = 0; iSamples < nSamples; ++iSamples)
			{
				if (iCurrentFoldTestSamples < nTestSamples && iSamples == samplesToBeTestOrder[iTestSamples])
				{
					matrixHeader[nTrainSamples + iCurrentFoldTestSamples] = "[TEST]_" + whole_cnt_matrix_output_header[iTestSamples];

					testSamplesText += whole_cnt_matrix_output_header[iSamples];
					if (iCurrentFoldTestSamples != nTestSamples - 1)
						testSamplesText += ", ";

					++iTestSamples;
					++iCurrentFoldTestSamples;
					continue;
				}
				matrixHeader[iSamples - iCurrentFoldTestSamples] = "[TRAIN]_" + whole_cnt_matrix_output_header[iSamples];
				trainSamplesText += whole_cnt_matrix_output_header[iSamples];
				if (iSamples != nSamples - 1)
					trainSamplesText += ", ";
			}

			auto logText = [&trainSamplesText, &testSamplesText](const std::string& alg, auto& top, auto& topCntMatrix, auto& topFasta)
			{
				Logger::Inst().Log("Info: generating " + alg + " cross-validation results to " + top + ".", 2);
				Logger::Inst().Log("Info: the file contains top k-mers with correlation values for " + trainSamplesText + " training, and " + testSamplesText + " testing sample(s).", 2);
				Logger::Inst().Log("Info: generating " + alg + " cross-validation results to " + topCntMatrix + ".", 2);
				Logger::Inst().Log("Info: the file contains counts matrix of top k-mers.", 2);
				Logger::Inst().Log("Info: generating " + alg + " cross-validation results to " + topFasta + ".", 2);
				Logger::Inst().Log("Info: the file contains top k-mers in FASTA format.", 2);
			};

			auto flush = [nSamples, iFold, matrixHeader, nCols, params, &tooLittleKmers, &nTopKmers, &nKmersExists, logText](StatisticsParams::CorrelationMethod correlationMethod, const std::string& header, const std::string& algorithmName, const auto& dataToFlush)
			{
				if (dataToFlush.empty())
					return;

				if (iFold == 0 && dataToFlush[iFold]->GetN() > dataToFlush[iFold]->Get().size())
				{
					assert((!tooLittleKmers && nKmersExists == 0) || nTopKmers == dataToFlush[iFold]->GetN());
					assert((!tooLittleKmers && nTopKmers == 0) || nKmersExists == dataToFlush[iFold]->Get().size());
					tooLittleKmers = true;

					nTopKmers = dataToFlush[iFold]->GetN();
					nKmersExists = dataToFlush[iFold]->Get().size();
				}

				const std::string top = params.statisticsParams.cvParams.getOuputFileNameTop(correlationMethod, nSamples, iFold, dataToFlush.size());
				const std::string topCntMatrix = params.statisticsParams.cvParams.getOuputFileNameTopCntMatrix(correlationMethod, nSamples, iFold, dataToFlush.size());
				const std::string topFasta = params.statisticsParams.cvParams.getOuputFileNameTopFasta(correlationMethod, nSamples, iFold, dataToFlush.size());
				flush_for(*dataToFlush[iFold],
					params.stage1Params.GetKmerLen(), nCols,
					top, { header },
					topCntMatrix, matrixHeader,
					topFasta);
				logText(algorithmName, top, topCntMatrix, topFasta);
			};

			flush(Pearson, "pearson", "Pearson", pearson);
			flush(Spearman, "spearman", "Spearman", spearman);
			flush(Kendall, "kendall", "Kendall Tau", kendall);
		}
		if (tooLittleKmers)
			Logger::Inst().Log("Warning: requested to store " + std::to_string(nTopKmers) + " CV k-mers with top p-values, but only " + std::to_string(nKmersExists) + " k-mers exist.");
	}

public:
	KeepNLargestCollectionCV(const std::vector<size_t>& samplesToBeTestOrder, size_t nTop, size_t nFolds, bool bPearson, bool bSpearman, bool bKendall) :
		samplesToBeTestOrder(&samplesToBeTestOrder)
	{
		if (nTop == 0)
			return;

		if (bPearson)
		{
			pearson.resize(nFolds);
			for (auto& it : pearson)
				it = std::make_unique<KeepTopNLargestABS_T>(nTop);
		}

		if (bSpearman)
		{
			spearman.resize(nFolds);
			for (auto& it : spearman)
				it = std::make_unique<KeepTopNLargestABS_T>(nTop);
		}

		if (bKendall)
		{
			kendall.resize(nFolds);
			for (auto& it : kendall)
				it = std::make_unique<KeepTopNLargestABS_T>(nTop);
		}
	}

	void addAllResultsForSingleStatistic(
		std::vector<std::unique_ptr<KeepTopNLargestABS_T>>& statisticResults,
		const std::string& kmerSeq,
		const kmcdb::CKmer<SIZE>& kmer,
		const std::vector<Statistics_T>& keys,
		const std::vector<VALUE_T>& allCounts) const
	{
		assert(statisticResults.size() == keys.size());
		const size_t nSamples = allCounts.size();
		const size_t nFolds = statisticResults.size();
		const size_t nTestSamples = nSamples / nFolds;
		const size_t nTrainSamples = nSamples - nTestSamples;

		std::vector<VALUE_T> counts(nSamples);

		size_t iTestSamples = 0;
		for (size_t iFold = 0; iFold < nFolds; ++iFold)
		{
			size_t iCurrentFoldTestSamples = 0;
			for (size_t iSamples = 0; iSamples < nSamples; ++iSamples)
			{
				if (iCurrentFoldTestSamples < nTestSamples && iSamples == (*samplesToBeTestOrder)[iTestSamples])
				{
					counts[nTrainSamples + iCurrentFoldTestSamples] = allCounts[iTestSamples];

					++iTestSamples;
					++iCurrentFoldTestSamples;
					continue;
				}
				counts[iSamples - iCurrentFoldTestSamples] = allCounts[iSamples];
			}

			statisticResults[iFold]->Add(Elem{ kmerSeq, kmer, keys[iFold], counts });
		}
	}

	void addPearson(
		const std::string& kmerSeq,
		const kmcdb::CKmer<SIZE>& kmer,
		const std::vector<Statistics_T>& keys,
		const std::vector<VALUE_T>& allCounts)
	{
		if (!pearson.empty())
			addAllResultsForSingleStatistic(pearson, kmerSeq, kmer, keys, allCounts);
	}

	void addSpearman(
		const std::string& kmerSeq,
		const kmcdb::CKmer<SIZE>& kmer,
		const std::vector<Statistics_T>& keys,
		const std::vector<VALUE_T>& allCounts)
	{
		if (!spearman.empty())
			addAllResultsForSingleStatistic(spearman, kmerSeq, kmer, keys, allCounts);
	}

	void addKendall(
		const std::string& kmerSeq,
		const kmcdb::CKmer<SIZE>& kmer,
		const std::vector<Statistics_T>& keys,
		const std::vector<VALUE_T>& allCounts)
	{
		if (!kendall.empty())
			addAllResultsForSingleStatistic(kendall, kmerSeq, kmer, keys, allCounts);
	}
};



template<unsigned SIZE, typename Statistics_T, typename VALUE_T, typename KeepNLargestCollection_T>
class KeepNLargestCollectionGlobal
{
private:
	std::mutex mtx;
	KeepNLargestCollection_T global;
public:
	void Add(KeepNLargestCollection_T& collection)
	{
		std::lock_guard lck(mtx);
		global.add(collection);
	}

	void Flush(const Params& params, const std::vector<std::string>& cnt_matrix_output_header)
	{
		global.flush(params, cnt_matrix_output_header);
	}
};



template<unsigned SIZE, typename Statistics_T, typename VALUE_T>
template<typename PRED>
void KeepNLargestCollectionBase<SIZE, Statistics_T, VALUE_T>::add_for(
	KeepNLargests<Elem, PRED>& src,
	KeepNLargests<Elem, PRED>& dest)
{
	assert(dest.GetN() == src.GetN()); //just to be sure that all source have the same N

	std::vector<Elem> data;
	src.Steal(data);
	for (auto& elem : data)
		dest.Add(std::move(elem));
}



template<unsigned SIZE, typename Statistics_T, typename VALUE_T>
template<typename PRED>
void KeepNLargestCollectionBase<SIZE, Statistics_T, VALUE_T>::flush_for(
	KeepNLargests<Elem, PRED>& to_flush,
	uint32_t first_col_len, size_t num_columns,
	const std::string& fname_top, const std::vector<std::string>& header_top,
	const std::string& fname_top_matrix, const std::vector<std::string>& header_top_matrix,
	const std::string& fname_top_fasta)
{
	TextFileWriter writer_top(fname_top, false);
	writer_top.StoreHeader(header_top);
	MatrixOutputBuffer<Statistics_T> buff_top(writer_top, first_col_len, 1);

	TextFileWriter writer_top_matrix(fname_top_matrix, false);
	writer_top_matrix.StoreHeader(header_top_matrix);
	MatrixOutputBuffer<VALUE_T> buff_top_matrix(writer_top_matrix, first_col_len, num_columns);

	TextFileWriter writer_fasta(fname_top_fasta, false);
	FastaOutputBuffer buff_top_fasta(writer_fasta, first_col_len);

	std::vector<Elem> data;
	to_flush.StealSorted(data, PRED{}); //could actually be Steal (no sorted), but lets keep it deterministic
	uint64_t outputKmerId = 0;
	bool containsNaNs = false;

	for (auto& elem : data)
	{
		if (std::isnan(elem.key))
			containsNaNs = true;

		buff_top.StoreKmer(elem.kmerSeq, elem.key);
		buff_top_matrix.StoreKmer(elem.kmerSeq, elem.counts);
		buff_top_fasta.StoreKmer(elem.kmerSeq, outputKmerId);
		++outputKmerId;
	}

	if (containsNaNs)
		Logger::Inst().Log("Warning: a file " + fname_top + " contains at least one NaN (not a number) value, thus results in " + fname_top_matrix + " and " + fname_top_fasta + " files may not be fully informative.");
}
