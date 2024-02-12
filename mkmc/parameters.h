#pragma once
#include <vector>
#include <string>
#include <thread>
#include <cstdint>
#include "kmc_core/kmc_runner.h"
#include "kmc_api/kmc_file.h"
#include "kmc_api/kmer_api.h"



enum class OutputFileType {FASTA, Matrix};

struct MKMCParams
{
	std::string tmpPath;
	std::vector<std::string> inputFiles;
	std::vector<std::string> kmcOutputFiles;
	std::vector<std::string> kmcTmpDirs;
	std::vector<std::string> toolsOutputFiles;
	OutputFileType outputFileType = OutputFileType::Matrix;

	uint32_t nThreads = std::thread::hardware_concurrency();
	uint32_t nKMCWorkers = 8;
	uint32_t maxRamGB = 16;
	std::string outputFile;
	bool dumpToFile = true;
	uint64_t dumpStepSize = static_cast<uint64_t>(1E9);
	const uint32_t count_symbols = 5;

	bool nKMCWorkersUserSet = false;

	bool keepTmpFiles = false;
};

struct FilterParams
{
	double minKmersAboveThresholdRatio = 0.0;
	uint32_t minCountThreshold = 1;
};

struct MutableParams
{
	bool tmpDirCreated = false;
};

struct Params
{
	KMC::Stage1Params stage1ParamsTemplate;
	KMC::Stage2Params stage2ParamsTemplate;
	MKMCParams mkmcParams;
	FilterParams filterParams;
	mutable MutableParams mutableParams;

	Params();
	void setKMCParams();
};