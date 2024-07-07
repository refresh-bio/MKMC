#pragma once

#include <cstdint>
#include <fstream>
#include <vector>
#include <string>
#include "KMCFileWrapper.h"
#include "kmc_dump/nc_utils.h"
#include "KmersSamplesStruct.h"
#include "parameters.h"
#include "kmcdb/bin_writers.h"


class MatrixFileGenerator
{
	kmcdb::BinWriterSortedPlain<uint64_t>* bin;
public:
	MatrixFileGenerator(const Params& params, uint32_t binId, kmcdb::BinWriterSortedPlain<uint64_t>* bin);

	template<typename KmersSamplesData_T>
	void writeKmer(const KmersSamplesData_T& kmersData);
};



class FASTAFileGenerator
{
	std::unique_ptr<char[]> str_kmer_buff;
	std::ofstream file;
	uint32_t k;

public:
	FASTAFileGenerator(const Params& params, uint32_t binId, kmcdb::BinWriterSortedPlain<uint64_t>* /*bin*/);

	template<typename KmersSamplesData_T>
	void writeKmer(const KmersSamplesData_T& kmersData);
};



template<typename Generator_T, typename... NextGenerators_T>
class PerformGenerate
{
	Generator_T generator;
	PerformGenerate<NextGenerators_T...> nextPerformGenerate;
public:
	PerformGenerate(const Params& params, uint32_t binId, kmcdb::BinWriterSortedPlain<uint64_t>* bin) :
		generator(params, binId, bin),
		nextPerformGenerate(params, binId, bin)
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
	PerformGenerate(const Params& params, uint32_t binId, kmcdb::BinWriterSortedPlain<uint64_t>* bin) :
		generator(params, binId, bin)
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
	bin->AddKmer(kmersData.minKmer, kmersData.kMersCounts.data());
}



template<typename KmersSamplesData_T>
void FASTAFileGenerator::writeKmer(const KmersSamplesData_T& kmersData)
{
	file << ">\n";
	kmersData.minKmer.to_string(k, str_kmer_buff.get());
	file << str_kmer_buff.get() << '\n';

	// k - number of written symbols
}
