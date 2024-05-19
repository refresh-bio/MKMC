#pragma once

#include <cstdint>
#include <fstream>
#include <vector>
#include <string>
#include "KMCFileWrapper.h"
#include "kmc_dump/nc_utils.h"
#include "KmersSamplesStruct.h"
#include "parameters.h"



class MatrixFileGenerator
{
	std::unique_ptr<char[]> str_kmer_buff;
	std::ofstream file;
	uint32_t k;

public:
	MatrixFileGenerator(const Params& params, uint32_t binId);

	template<typename KmersSamplesData_T>
	void writeKmer(const KmersSamplesData_T& kmersData);
};



class FASTAFileGenerator
{
	std::unique_ptr<char[]> str_kmer_buff;
	std::ofstream file;
	uint32_t k;

public:
	FASTAFileGenerator(const Params& params, uint32_t binId);

	template<typename KmersSamplesData_T>
	void writeKmer(const KmersSamplesData_T& kmersData);
};



template<typename Generator_T, typename... NextGenerators_T>
class PerformGenerate
{
	Generator_T generator;
	PerformGenerate<NextGenerators_T...> nextPerformGenerate;
public:
	PerformGenerate(const Params& params, uint32_t binId) :
		generator(params, binId),
		nextPerformGenerate(params, binId)
	{}

	template<typename KmersSamplesData_T>
	void writeKmer(const KmersSamplesData_T& kmersData)
	{
		generator.writeKmer(kmersData);
		nextPerformGenerate.writeKmer(kmersData);
	}
};



template<typename Generator_T>
class PerformGenerate<Generator_T>
{
	Generator_T generator;
public:
	PerformGenerate(const Params& params, uint32_t binId) :
		generator(params, binId)
	{}

	template<typename KmersSamplesData_T>
	void writeKmer(const KmersSamplesData_T& kmersData)
	{
		generator.writeKmer(kmersData);
	}
};



template<typename KmersSamplesData_T>
void MatrixFileGenerator::writeKmer(const KmersSamplesData_T& kmersData)
{
	kmersData.minKmer.to_string(k, str_kmer_buff.get());
	uint32_t pos = k;
	for (uint64_t count : kmersData.kMersCounts)
	{
		str_kmer_buff[pos++] = '\t';
		uint32_t shift = CNumericConversions::Int2PChar(count, reinterpret_cast<uchar*>(str_kmer_buff.get()) + pos);
		pos += shift;
	}

	str_kmer_buff[pos] = '\n';
	str_kmer_buff[pos + 1] = '\0';
	file << str_kmer_buff.get();

	// pos - number of written symbols
}



template<typename KmersSamplesData_T>
void FASTAFileGenerator::writeKmer(const KmersSamplesData_T& kmersData)
{
	file << ">\n";
	kmersData.minKmer.to_string(k, str_kmer_buff.get());
	file << str_kmer_buff.get() << '\n';

	// k - number of written symbols
}
