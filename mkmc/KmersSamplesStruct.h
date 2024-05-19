#pragma once

#include <vector>



template <unsigned _SIZE>
struct KmersSamplesStruct
{
	static const unsigned SIZE = _SIZE;

	CKmer<_SIZE> minKmer;
	const std::vector<uint64_t>& kMersCounts;
};
