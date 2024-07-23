#pragma once

#include <vector>



template <unsigned _SIZE>
struct KmersSamplesStruct
{
	static const unsigned SIZE = _SIZE;

	kmcdb::CKmer<_SIZE> kmer;
	const std::vector<uint64_t>& kMersCounts;
};
