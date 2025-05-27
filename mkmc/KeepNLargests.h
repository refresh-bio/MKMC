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
			return std::make_pair(std::abs(lhs.key), rhs.kmer) > std::make_pair(std::abs(rhs.key), lhs.kmer);
		}
	};

	struct Greater
	{
		bool operator()(const KeepTopElem& lhs, const KeepTopElem& rhs)
		{
			return std::make_pair(lhs.key, rhs.kmer) > std::make_pair(rhs.key, lhs.kmer);
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
	friend struct KeepNLargestCollectionGlobal;

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
	void add_impl(std::unique_ptr<KeepNLargests<Elem, PRED>>& src, std::unique_ptr<KeepNLargests<Elem, PRED>>& dest)
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
		flush_for(*pearson,
			params.stage1Params.GetKmerLen(), params.mkmcParams.samples.size(),
			params.mkmcParams.outputFilePearsonTop, { "pearson" },
			params.mkmcParams.outputFilePearsonTopCntMatrix, cnt_matrix_output_header,
			params.mkmcParams.outputFilePearsonTopFasta);

		flush_for(*spearman,
			params.stage1Params.GetKmerLen(), params.mkmcParams.samples.size(),
			params.mkmcParams.outputFileSpearmanTop, { "spearman" },
			params.mkmcParams.outputFileSpearmanTopCntMatrix, cnt_matrix_output_header,
			params.mkmcParams.outputFileSpearmanTopFasta);

		flush_for(*kendall,
			params.stage1Params.GetKmerLen(), params.mkmcParams.samples.size(),
			params.mkmcParams.outputFileKendallTop, { "kendall" },
			params.mkmcParams.outputFileKendallTopCntMatrix, cnt_matrix_output_header,
			params.mkmcParams.outputFileKendallTopFasta);

		flush_for(*entropy,
			params.stage1Params.GetKmerLen(), params.mkmcParams.samples.size(),
			params.mkmcParams.outputFileEntropyTop, { "entropy" },
			params.mkmcParams.outputFileEntropyTopCntMatrix, cnt_matrix_output_header,
			params.mkmcParams.outputFileEntropyTopFasta);

		flush_for(*snr,
			params.stage1Params.GetKmerLen(), params.mkmcParams.samples.size(),
			params.mkmcParams.outputFileSNRTop, { "snr" },
			params.mkmcParams.outputFileSNRTopCntMatrix, cnt_matrix_output_header,
			params.mkmcParams.outputFileSNRTopFasta);

		flush_for(*unnormalizedSnr,
			params.stage1Params.GetKmerLen(), params.mkmcParams.samples.size(),
			params.mkmcParams.outputFileUnnormalizedSNRTop, { "snr_for_unnormalized" },
			params.mkmcParams.outputFileUnnormalizedSNRTopCntMatrix, cnt_matrix_output_header,
			params.mkmcParams.outputFileUnnormalizedSNRTopFasta);

		flush_for(*dids,
			params.stage1Params.GetKmerLen(), params.mkmcParams.samples.size(),
			params.mkmcParams.outputFileDIDSTop, { "dids" },
			params.mkmcParams.outputFileDIDSTopCntMatrix, cnt_matrix_output_header,
			params.mkmcParams.outputFileDIDSTopFasta);
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
	for (auto& elem : data)
	{
		buff_top.StoreKmer(elem.kmerSeq, elem.key);
		buff_top_matrix.StoreKmer(elem.kmerSeq, elem.counts);
		buff_top_fasta.StoreKmer(elem.kmerSeq, outputKmerId);
		++outputKmerId;
	}
}
