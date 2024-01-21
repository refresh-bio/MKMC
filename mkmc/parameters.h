#pragma once
#include <vector>
#include <string>
#include <thread>
#include "kmc_core/kmc_runner.h"
#include "kmc_api/kmc_file.h"
#include "kmc_api/kmer_api.h"

struct MKMCParams
{
	std::vector<std::string> inputFiles;
	std::vector<std::string> kmcOutputFiles;
	std::vector<std::string> toolsOutputFiles;
	uint32_t nThreads = std::thread::hardware_concurrency();
	uint32_t nKMCWorkers = 2;
	uint32_t maxRamGB = 12;
	std::string outputFile;
	bool dumpToFile = false;
	uint64_t dumpStepSize = static_cast<uint64_t>(1E9);
	double minKmersPresenceThreshold = 0.0;
	const uint32_t count_symbols = 5;
};

struct Params
{
	KMC::Stage1Params stage1ParamsTemplate;
	KMC::Stage2Params stage2ParamsTemplate;
	MKMCParams mkmcParams;

	Params();
	void setKMCParams();
};