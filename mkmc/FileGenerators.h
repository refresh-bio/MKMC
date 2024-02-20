#pragma once

#include <cstdint>
#include <fstream>
#include <vector>
#include <string>
#include "KMCFileWrapper.h"



class MatrixFileGenerator
{
	std::unique_ptr<char[]> str_kmer_buff;
	std::ostream& file;
	uint32_t k;

public:
	MatrixFileGenerator(std::ostream& file, const std::vector<std::string>& samples, uint32_t countSymbols, uint32_t k);

	uint32_t writeKmer(KMCFileWrapper::kmer_t& kmer, const std::vector<size_t>& kMersCounts);
};



class FASTAFileGenerator
{
	std::unique_ptr<char[]> str_kmer_buff;
	std::ostream& file;
	uint32_t k;

public:
	FASTAFileGenerator(std::ostream& file, const std::vector<std::string>& samples, uint32_t countSymbols, uint32_t k);

	uint32_t writeKmer(KMCFileWrapper::kmer_t& kmer, const std::vector<size_t>& kMersCounts);
};
