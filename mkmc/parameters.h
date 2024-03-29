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
	KMC::InputFileType inputFileType = KMC::InputFileType::FASTQ;
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

struct FilterParams
{
	double minKmersAboveThresholdRatio = 0.0;
	uint32_t minCountThreshold = 1;

	bool filterKmersSequences = false;
	std::string inputKmersSequencesToFilterOut; // input file name
	std::string kmersSequencesToFilterOutDB;    // KMC file name
};

struct MutableParams
{
	bool tmpDirCreated = false;

	mutable bool createdFastaFile = false;
	mutable std::string kmersSequencesToFilterOut; // file containing k-mers to count (input or generated)
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