#pragma once

//#include "refresh/conversions/libs/conversions.h" - has to be included by kmcdb
#include <cstdint>
#include <fstream>
#include <vector>
#include <string>
#include "KMCFileWrapper.h"
#include "kmc_dump/nc_utils.h"
#include "KmersSamplesStruct.h"
#include "parameters.h"
#include "Logger.h"
#include "kmcdb/kmcdb.h"
#include "kmcdb/bin_writers.h"
#include "TextFileWritingUtilities.h"



/* Generators are designed to be assigned to threads. Theoretically may be also assigned to single tasks (bins).
 * closeWriter functions close writers manually to prevent waiting for closing to the end of the program,
 * as writers are static (and common for all threads).
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
	void writeKmer(const KmersSamplesData_T& kmersData, uint64_t kmerIdInBin);

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
					params.mkmcParams.outputMatrixBinFile,
					"",
					sample_names);
			}
			catch (const std::exception& ex)
			{
				Logger::Inst().Log(std::string("Error: ") + ex.what());
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
	static TextFileWriter* dumpWriter;
	static bool writerWasOpened;

	MatrixOutputBuffer<uint64_t> outputBuffer;
	const uint32_t kmerLength;
	std::string kmerSeqBuf;

public:
	MatrixFileGenerator(const Params& params) :
		outputBuffer(*dumpWriter, kmerLength, params.mkmcParams.samples.size()),
		kmerLength(params.stage1Params.GetKmerLen()),
		kmerSeqBuf(params.stage1Params.GetKmerLen(), ' ')
	{}

	void setBinId(uint32_t binId){}

	template<typename KmersSamplesData_T>
	void writeKmer(const KmersSamplesData_T& kmersData, uint64_t kmerIdInBin);

	static void initWriter(const Params& params, const kmcdb::Config& config)
	{
		assert(!writerWasOpened);
		if (dumpWriter == nullptr)
		{
			std::vector<std::string> sampleNames;
			sampleNames.reserve(params.mkmcParams.samples.size());
			for (const auto& sample : params.mkmcParams.samples)
				sampleNames.push_back(sample.name);

			dumpWriter = new TextFileWriter(params.mkmcParams.outputMatrixFile, (params.mkmcParams.nThreads > 1));
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
	static TextFileWriter* dumpWriter;
	static bool writerWasOpened;

	uint32_t binId;

	FastaOutputBuffer outputBuffer;
	uint64_t kmerLength;
	std::string kmerSeqBuf;
public:
	FASTAFileGenerator(const Params& params) :
		binId(std::numeric_limits<uint32_t>::max()),
		outputBuffer(*dumpWriter, params.stage1Params.GetKmerLen()),
		kmerLength(params.stage1Params.GetKmerLen()),
		kmerSeqBuf(params.stage1Params.GetKmerLen(), ' ')
	{}

	void setBinId(uint32_t binId)
	{
		this->binId = binId;
	}

	template<typename KmersSamplesData_T>
	void writeKmer(const KmersSamplesData_T& kmersData, uint64_t kmerIdInBin);

	static void initWriter(const Params& params, const kmcdb::Config& config)
	{
		assert(!writerWasOpened);
		if (dumpWriter == nullptr)
		{
			dumpWriter = new TextFileWriter(params.mkmcParams.outputFASTAFile, (params.mkmcParams.nThreads > 1));
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
	void writeKmer(const KmersSamplesData_T& kmersData, uint64_t kmerIdInBin)
	{
		generator.writeKmer(kmersData, kmerIdInBin);
		nextPerformGenerate.writeKmer(kmersData, kmerIdInBin);
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
	void writeKmer(const KmersSamplesData_T& kmersData, uint64_t kmerIdInBin)
	{
		generator.writeKmer(kmersData, kmerIdInBin);
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
void BinFileGenerator::writeKmer(const KmersSamplesData_T& kmersData, uint64_t kmerIdInBin)
{
	assert(kmcBinDBWriter != nullptr);
	kmcBinDBWriter->AddKmer(kmersData.kmer, kmersData.kMersCounts.data());
}



template<typename KmersSamplesData_T>
void MatrixFileGenerator::writeKmer(const KmersSamplesData_T& kmersData, uint64_t kmerIdInBin)
{
	kmersData.kmer.to_string(kmerLength, kmerSeqBuf.data());
	outputBuffer.StoreKmer(kmerSeqBuf, kmersData.kMersCounts);
}



template<typename KmersSamplesData_T>
void FASTAFileGenerator::writeKmer(const KmersSamplesData_T& kmersData, uint64_t kmerIdInBin)
{
	kmersData.kmer.to_string(kmerLength, kmerSeqBuf.data());
	outputBuffer.StoreKmer(kmerSeqBuf, binId, kmerIdInBin);
}
