#pragma once
#include <vector>
#include <string>
#include <thread>
#include <cstdint>
#include "kmc_core/kmc_runner.h"
#include "kmc_api/kmc_file.h"
#include "kmc_api/kmer_api.h"
#undef small
#include "lib/statistics/lib/statistics_normalization.h"



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
	std::vector<std::string> outputFASTAFiles;
	std::vector<std::string> outputMatrixFiles;
	std::vector<OutputFileType> outputFileTypes = { OutputFileType::Matrix };


	std::vector<std::string> outputFilesNormFrequency;
	std::vector<std::string> outputFilesNormQuantile;

	std::vector<std::string> outputFilesPearson;
	std::vector<std::string> outputFilesSpearman;
	std::vector<std::string> outputFilesKendall;

	uint32_t nThreads = (std::min)(16U, std::thread::hardware_concurrency());
	uint32_t nKMCWorkers = 4;
	uint32_t maxRamGB = 16;
	const uint32_t countSymbols = 10; // for 4G

	const uint32_t sigToBinMapStatsPercentage = 5;
	uint32_t nKMCBins = 512;

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

struct StatisticsParams
{
	using NormalizationMethod = refresh::normalization_base<size_t, double>::method_t;
	using NormalizationLearning = refresh::normalization_learn<size_t, double>;

	bool generateStatistics = false;
	std::string phenotypeFile;

	std::string normFrequencyFileTmp = "frequencyDump";
	std::string normQuantileFileTmp = "quantileDump";
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
	StatisticsParams statisticsParams;

	Params();
	void setKMCParams();
};