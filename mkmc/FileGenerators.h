#pragma once

//#include "lib/refresh/conversions/libs/conversions.h" - has to be included by kmcdb
#include <cstdint>
#include <fstream>
#include <vector>
#include <string>
#include "KMCFileWrapper.h"
#include "kmc_dump/nc_utils.h"
#include "KmersSamplesStruct.h"
#include "parameters.h"
#include "kmcdb/kmcdb.h"
#include "kmcdb/bin_writers.h"
#include "DumpWriter.h"



/* Generators are designed to be assigned to threads. Theoretically may be also assigned to single tasks (bins).
 *  closeWriter functions close writers manually to prevent waiting for closing to the end of the program.
 */
class BinFileGenerator
{
	static kmcdb::WriterSortedPlain<uint64_t>* kmcDBWriter;
	static bool writerWasOpened;

	kmcdb::BinWriterSortedPlain<uint64_t>* kmcBinDBWriter;
	uint64_t kmerLength;
public:
	BinFileGenerator(const Params& params) :
		kmcBinDBWriter(nullptr),
		kmerLength(params.stage1Params.GetKmerLen())
	{}

	void setBinId(uint32_t binId)
	{
		kmcBinDBWriter = kmcDBWriter->GetBin(binId);
	}

	template<typename KmersSamplesData_T>
	void writeKmer(const KmersSamplesData_T& kmersData);

	static void initWriter(const Params& params, const kmcdb::Config& config)
	{
		assert(!writerWasOpened);
		if (kmcDBWriter == nullptr)
		{
			std::vector<std::string> sample_names{};
			sample_names.reserve(params.mkmcParams.samples.size());

			for (const auto& sample : params.mkmcParams.samples)
				sample_names.push_back(sample.name);

			kmcdb::ConfigSortedPlain representation_config{};

			try
			{
				kmcDBWriter = new kmcdb::WriterSortedPlain<uint64_t>(
					config,
					representation_config,
					params.mkmcParams.outputBinFile,
					"",
					sample_names);
			}
			catch (const std::exception& ex)
			{
				std::cerr << "Error: " << ex.what() << "\n";
				exit(1);
			}
			writerWasOpened = true;
		}
	}

	static void closeWriter()
	{
		delete kmcDBWriter;
		kmcDBWriter = nullptr;
	}
};



class MatrixFileGenerator
{
	static DumpWriter* dumpWriter;
	static bool writerWasOpened;

	OutputBuffer outputBuffer;
	uint32_t kmerLength;

	size_t getMaxLineLength(const Params& params) const
	{
		return params.stage1Params.GetKmerLen() + 1 + params.mkmcParams.samples.size() * (refresh::numeric_conversion_max_length<decltype(params.stage1Params.GetKmerLen())>() + 1);
	}
public:
	MatrixFileGenerator(const Params& params) :
		outputBuffer(*dumpWriter, getMaxLineLength(params)),
		kmerLength(params.stage1Params.GetKmerLen())
	{}

	void setBinId(uint32_t binId){}

	template<typename KmersSamplesData_T>
	void writeKmer(const KmersSamplesData_T& kmersData);

	static void initWriter(const Params& params, const kmcdb::Config& config)
	{
		assert(!writerWasOpened);
		if (dumpWriter == nullptr)
		{
			std::vector<std::string> sampleNames;
			sampleNames.reserve(params.mkmcParams.samples.size());
			for (const auto& sample : params.mkmcParams.samples)
				sampleNames.push_back(sample.name);

			dumpWriter = new DumpWriter(params.mkmcParams.outputMatrixFile, (params.mkmcParams.nThreads > 1));
			dumpWriter->StoreHeader(sampleNames);
			writerWasOpened = true;
		}
	}
	static void closeWriter()
	{
		delete dumpWriter;
		dumpWriter = nullptr;
	}
};



class FASTAFileGenerator
{
	static DumpWriter* dumpWriter;
	static bool writerWasOpened;

	OutputBuffer outputBuffer;
	uint64_t kmerLength;
public:
	FASTAFileGenerator(const Params& params) :
		outputBuffer(*dumpWriter, params.stage1Params.GetKmerLen()),
		kmerLength(params.stage1Params.GetKmerLen())
	{}

	void setBinId(uint32_t binId) {}

	template<typename KmersSamplesData_T>
	void writeKmer(const KmersSamplesData_T& kmersData);

	static void initWriter(const Params& params, const kmcdb::Config& config)
	{
		assert(!writerWasOpened);
		if (dumpWriter == nullptr)
		{
			dumpWriter = new DumpWriter(params.mkmcParams.outputFASTAFile, (params.mkmcParams.nThreads > 1));
			writerWasOpened = true;
		}
	}
	static void closeWriter()
	{
		delete dumpWriter;
		dumpWriter = nullptr;
	}
};



template<typename Generator_T, typename... NextGenerators_T>
class PerformGenerate
{
	Generator_T generator;
	PerformGenerate<NextGenerators_T...> nextPerformGenerate;
public:
	PerformGenerate(const Params& params) :
		generator(params),
		nextPerformGenerate(params)
	{}

	void setBinId(uint32_t binId)
	{
		generator.setBinId(binId);
		nextPerformGenerate.setBinId(binId);
	}

	template<typename KmersSamplesData_T>
	void writeKmer(const KmersSamplesData_T& kmersData)
	{
		generator.writeKmer(kmersData);
		nextPerformGenerate.writeKmer(kmersData);
	}

	static void initWriters(const Params& params, const kmcdb::Config& config)
	{
		PerformGenerate<NextGenerators_T...>::initWriters(params, config);
		Generator_T::initWriter(params, config);
	}

	static void closeWriters()
	{
		PerformGenerate<NextGenerators_T...>::closeWriters();
		Generator_T::closeWriter();
	}
};



template<typename Generator_T>
class PerformGenerate<Generator_T>
{
	Generator_T generator;
public:
	PerformGenerate(const Params& params) :
		generator(params)
	{}

	void setBinId(uint32_t binId)
	{
		generator.setBinId(binId);
	}

	template<typename KmersSamplesData_T>
	void writeKmer(const KmersSamplesData_T& kmersData)
	{
		generator.writeKmer(kmersData);
	}

	static void initWriters(const Params& params, const kmcdb::Config& config)
	{
		Generator_T::initWriter(params, config);
	}

	static void closeWriters()
	{
		Generator_T::closeWriter();
	}
};



template<typename KmersSamplesData_T>
void BinFileGenerator::writeKmer(const KmersSamplesData_T& kmersData)
{
	assert(kmcBinDBWriter != nullptr);
	kmcBinDBWriter->AddKmer(kmersData.kmer, kmersData.kMersCounts.data());
}



template<typename KmersSamplesData_T>
void MatrixFileGenerator::writeKmer(const KmersSamplesData_T& kmersData)
{
	auto storeMethod = []<unsigned SIZE, typename VALUE_T>(const kmcdb::CKmer<SIZE>& kmer, uint64_t kmer_len, const std::vector<VALUE_T>& cnts, char* out) -> size_t
	{
		kmer.to_string(kmer_len, out, '\t');
		size_t res = kmer_len + 1;
		out += kmer_len + 1;

		auto store_single_value = [&](const VALUE_T& val, char term)
		{
			size_t r{};
			if constexpr (std::is_integral_v<VALUE_T>)
				r = refresh::int_to_pchar(val, out, term);
			else if constexpr (std::is_floating_point_v<VALUE_T>)
				r = refresh::real_to_pchar(val, out, 6, term);
			else
			{
				static_assert(!sizeof(VALUE_T), "Unsupported type");
			}
			out += r;
			res += r;
		};

		for (size_t i = 0; i < cnts.size() - 1; ++i)
			store_single_value(cnts[i], '\t');
		store_single_value(cnts.back(), '\n');

		return res;
	};

	outputBuffer.StoreKmer(kmersData.kmer, kmerLength, kmersData.kMersCounts, storeMethod);
}



template<typename KmersSamplesData_T>
void FASTAFileGenerator::writeKmer(const KmersSamplesData_T& kmersData)
{
	auto storeMethod = []<unsigned SIZE, typename VALUE_T>(const kmcdb::CKmer<SIZE>& kmer, uint64_t kmer_len, const std::vector<VALUE_T>& cnts, char* out) -> size_t
	{
		out[0] = '>';
		out[1] = '\n';
		size_t res = 2;
		out += 2;

		kmer.to_string(kmer_len, out, '\n');
		res += kmer_len + 1;
		out += kmer_len + 1;

		return res;
	};

	outputBuffer.StoreKmer(kmersData.kmer, kmerLength, kmersData.kMersCounts, storeMethod);
}
