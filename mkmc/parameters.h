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
	std::vector<std::string> samples;
	std::vector<std::vector<std::string>> inputFilesPerSample;
	std::vector<std::string> kmcOutputFiles;
	std::vector<std::string> kmcTmpDirs;

	std::string outputFilesTemplate;
	std::vector<std::string> outputFiles;
	OutputFileType outputFileType = OutputFileType::Matrix;

	uint32_t nThreads = std::thread::hardware_concurrency();
	uint32_t nDumpThreads = (std::min)(16U, std::thread::hardware_concurrency());
	uint32_t nKMCWorkers = 8;
	uint32_t maxRamGB = 16;
	uint64_t dumpStepSize = static_cast<uint64_t>(1E9);
	bool splitDumpOutput = false;
	const uint32_t count_symbols = 10; // for 4G

	const uint32_t sigToBinMapStatsPercentage = 5;

	bool nKMCWorkersUserSet = false;

	bool keepTmpFiles = false;

	int verbosity_level = 0;
};

struct KMCToolsParams
{
	uint32_t nThreads = 1;
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
	MKMCParams mkmcParams;
	mutable MutableParams mutableParams;

	KMC::Stage1Params stage1Params;
	KMC::Stage2Params stage2Params;

	FilterParams filterParams;

	Params();
	void setKMCParams();
};