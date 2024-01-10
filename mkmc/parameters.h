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
	std::vector<std::string> tmpFiles;
	uint32_t nThreads = std::thread::hardware_concurrency();
	uint32_t nKMCWorkers = 2;
	uint32_t maxRamGB = 12;
	std::string outputFile;
};

struct Params
{
	KMC::Stage1Params stage1ParamsTemplate;
	KMC::Stage2Params stage2ParamsTemplate;
	MKMCParams mkmcParams;

	Params();
	void setKMCParams();
};