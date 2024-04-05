#pragma once

#include <cstdint>
#include <fstream>
#include <vector>
#include <string>
#include "KMCFileWrapper.h"
#include "../kmc/kmc_dump/nc_utils.h"
#include "KmersSamplesStruct.h"



class MatrixFileGenerator
{
	std::unique_ptr<char[]> str_kmer_buff;
	std::ostream& file;
	uint32_t k;

public:
	MatrixFileGenerator(std::ostream& file, const std::vector<std::string>& samples, uint32_t countSymbols, uint32_t k);

	template<typename KmersSamplesData_T>
	uint32_t writeKmer(const KmersSamplesData_T& kmersData);
};



class FASTAFileGenerator
{
	std::unique_ptr<char[]> str_kmer_buff;
	std::ostream& file;
	uint32_t k;

public:
	FASTAFileGenerator(std::ostream& file, const std::vector<std::string>& samples, uint32_t countSymbols, uint32_t k);

	template<typename KmersSamplesData_T>
	uint32_t writeKmer(const KmersSamplesData_T& kmersData);
};



template<typename KmersSamplesData_T>
uint32_t MatrixFileGenerator::writeKmer(const KmersSamplesData_T& kmersData)
{
	kmersData.minKmer.to_string(k, str_kmer_buff.get());
	uint32_t pos = k;
	for (size_t count : kmersData.kMersCounts)
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


template<typename KmersSamplesData_T>
uint32_t FASTAFileGenerator::writeKmer(const KmersSamplesData_T& kmersData)
{
	file << ">\n";
	kmersData.minKmer.to_string(k, str_kmer_buff.get());
	file << str_kmer_buff.get() << '\n';

	return k;
}
