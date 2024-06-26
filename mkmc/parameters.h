#pragma once
#include <vector>
#include <string>
#include <thread>
#include <cstdint>
#include "kmc_core/kmc_runner.h"
#include "kmc_api/kmc_file.h"
#include "kmc_api/kmer_api.h"
#undef small
#include "refresh/statistics/lib/statistics_normalization.h"



enum class OutputFileType {FASTA, Matrix};

struct Sample
{
	std::string name;
	std::vector<std::string> inputFiles;
};

struct MKMCParams
{
	std::string inputFileName;
	std::string tmpPath;
	std::vector<Sample> samples;

	KMC::InputFileType inputFileType = KMC::InputFileType::FASTQ; // also default "fq" value in input parameters
	std::vector<std::string> kmcOutputFiles;
	std::vector<std::string> kmcTmpDirs;

	std::string outputFilesTemplate;
	std::vector<std::string> outputFASTAFiles;
	std::vector<std::string> outputMatrixFiles; //mkokot_TODO: this will be probably to be removed
	std::vector<OutputFileType> outputFileTypes = { OutputFileType::Matrix }; // also default "matrix" value in input parameters


	std::vector<std::string> outputFilesNorm;

	std::vector<std::string> outputFilesPearson;
	std::vector<std::string> outputFilesSpearman;
	std::vector<std::string> outputFilesKendall;

	uint32_t nThreads = (std::min)(16U, std::thread::hardware_concurrency());
	uint32_t nKMCWorkers = 4;
	uint32_t maxRamGB = 16;
	bool maxRamGBUserDefined = false;
	const uint32_t countSymbols = 10; // for 4G

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
	using NormalizationMethod = refresh::normalization_base<uint64_t, double>::method_t;
	using NormalizationLearning = refresh::normalization_learn<uint64_t, double>;

	enum class CorrelationMethod { Pearson, Spearman, Kendall };

	bool generateNormalization = false;
	NormalizationMethod normalizationMethod;

	std::string phenotypeFile;
	std::vector<CorrelationMethod> correlationMethods;

	std::string normFrequencyFileTmp = "frequencyDump";
	std::string normQuantileFileTmp = "quantileDump";

	std::string statsNOutputKmers = "nKmers";
};

struct MutableParams
{
	bool tmpDirCreated = false;

	mutable bool createdFastaFile = false;
	mutable std::string kmersSequencesToFilterOut; // file containing k-mers to count (input or generated)
};

struct DefaultKMCParams
{
	const uint32_t k = 25;
	const uint64_t ci = 1;
	const uint64_t cx = static_cast<uint32_t>(4e9);
	const uint64_t cs = 65535;
};

struct Params
{
	MKMCParams mkmcParams;
	mutable MutableParams mutableParams;

	DefaultKMCParams defaultKMCParams;

	KMC::Stage1Params stage1Params;
	KMC::Stage2Params stage2Params;

	FilterParams filterParams;
	StatisticsParams statisticsParams;

	Params();
	void generateTempAndOutputFilesNames();
	bool readAdditionalParamsFromFiles();
	void adjustKMCPerformanceParams();
	void adjustAnotherParams();
};