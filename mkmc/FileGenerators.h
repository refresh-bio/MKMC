#pragma once

#include <cstdint>
#include <fstream>
#include <vector>
#include <string>
#include "KMCFileWrapper.h"
#include "../kmc/kmc_dump/nc_utils.h"

class MatrixFileGenerator
{
	std::unique_ptr<char[]> str_kmer_buff;
	std::ostream& file;
	uint32_t k;

public:
	MatrixFileGenerator(std::ostream& file, const std::vector<std::string>& samples, uint32_t countSymbols, uint32_t k);

	template<unsigned SIZE>
	uint32_t writeKmer(CKmer<SIZE>& kmer, const std::vector<size_t>& kMersCounts);
};



class FASTAFileGenerator
{
	std::unique_ptr<char[]> str_kmer_buff;
	std::ostream& file;
	uint32_t k;

public:
	FASTAFileGenerator(std::ostream& file, const std::vector<std::string>& samples, uint32_t countSymbols, uint32_t k);

	template<unsigned SIZE>
	uint32_t writeKmer(CKmer<SIZE>& kmer, const std::vector<size_t>& kMersCounts);
};

template<unsigned SIZE>
uint32_t MatrixFileGenerator::writeKmer(CKmer<SIZE>& kmer, const std::vector<size_t>& kMersCounts)
{
	kmer.to_string(k, str_kmer_buff.get());
	uint32_t pos = k;
	for (size_t count : kMersCounts)
	{
		str_kmer_buff[pos++] = '\t';
		uint32_t shift = CNumericConversions::Int2PChar(count, reinterpret_cast<uchar*>(str_kmer_buff.get()) + pos);
		pos += shift;
	}

	str_kmer_buff[pos] = '\n';
	str_kmer_buff[pos + 1] = '\0';
	file << str_kmer_buff.get();

	return pos;
}


template<unsigned SIZE>
uint32_t FASTAFileGenerator::writeKmer(CKmer<SIZE>& kmer, const std::vector<size_t>& kMersCounts)
{
	file << ">\n";
	kmer.to_string(k, str_kmer_buff.get());
	file << str_kmer_buff.get() << '\n';

	return k;
}
