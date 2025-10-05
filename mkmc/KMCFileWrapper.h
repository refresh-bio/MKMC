#pragma once

#include <string>
#include <cstdint>
#include <iostream>
#include "kmcdb/bin_readers.h"


//mkokot_TODO: rename this class
template<unsigned SIZE>
class KMCFileWrapper
{
public:

private:
	kmcdb::BinReaderSortedWithLUTForListing<uint64_t>* bin;
	size_t tot_kmers;
	size_t cur_kmer_no{};
	kmcdb::CKmer<SIZE> cur;
	uint64_t cur_count;
public:
	KMCFileWrapper(KMCFileWrapper&&) = default;
	KMCFileWrapper& operator=(KMCFileWrapper&&) = default;
	KMCFileWrapper(kmcdb::BinReaderSortedWithLUTForListing<uint64_t>* bin);

	bool Finished() const
	{
		return cur_kmer_no > tot_kmers;
	}
	const kmcdb::CKmer<SIZE>& First() const
	{
		return cur;
	}
	uint64_t FirstCount() const
	{
		return cur_count;
	}
	size_t GetTotKmers() const
	{
		return tot_kmers;
	}
	void Next();
};

template<unsigned SIZE>
KMCFileWrapper<SIZE>::KMCFileWrapper(kmcdb::BinReaderSortedWithLUTForListing<uint64_t>* bin) :
	bin(bin)
{
	tot_kmers = bin->GetBinMetadata().total_kmers;
	cur.clear();

	Next();
}


template<unsigned SIZE>
inline void KMCFileWrapper<SIZE>::Next()
{
	if (!bin->NextKmer(cur, &cur_count))
	{
		if (cur_kmer_no != tot_kmers)
		{
			std::cerr << "Error: critical, this should not happen, details: " << __FILE__ << "(" << __LINE__ << ")." << std::endl;
			exit(1);
		}
	}
	++cur_kmer_no;
}
