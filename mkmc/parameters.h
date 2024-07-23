#pragma once
#include <vector>
#include <string>
#include <thread>
#include <cstdint>
#include "PhenotypeReaders.h"
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
	std::string outputBinFile;
	std::string outputStatsBinFile;
	std::string outputFASTAFile;
	std::string outputMatrixFile;
	std::vector<OutputFileType> outputFileTypes;

	std::string normStatsBinFile;

	std::string outputFileNorm;

	std::string outputFilePearson;
	std::string outputFileSpearman;
	std::string outputFileKendall;

	std::string outputFileEntropy;

	std::string outputFileTTest;
	std::string outputFileTTestCor;
	std::string outputFileSNR;
	std::string outputFileWilcoxonRankSum;
	std::string outputFileWilcoxonRankSumCor;
	std::string outputFileDIDS;
	std::string outputFileANOVA;
	std::string outputFileANOVACor;

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

	bool generateNormalization = false;
	NormalizationMethod normalizationMethod;

	enum class CorrelationMethod { Pearson, Spearman, Kendall };
	std::vector<CorrelationMethod> correlationMethods;

	enum class DifferentialAnalysisMethod { TTest, SNR, WilcoxonRankSum, DIDS, ANOVA };
	std::vector<DifferentialAnalysisMethod> classificationMethods;

	bool correctPvalues = false;
	enum class DifferentialAnalysisCorrectionMethod { Bonferroni, HolmBonferroni, BenjaminiHochberg, BenjaminiYekutieli };
	DifferentialAnalysisCorrectionMethod classificationPValueCorrection;

	bool generateEntropy = false;
	inline const static std::string normDeseq2StreamName = "norm_deseq2";
	inline const static std::string normFrequencyStreamName = "norm_frequency";
	inline const static std::string normQuantileStreamName = "norm_quantile";
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

struct Phenotypes
{
	PhenotypeReader<int64_t> correlationPhenotype;
	DifferentialAnalysisPhenotypeReader differentialAnalysisPhenotype; 
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

	Phenotypes phenotypes;

	Params();
	void generateTempAndOutputFilesNames();
	bool readAdditionalParamsFromFiles();
	void adjustKMCPerformanceParams();
	void adjustAnotherParams();
	void readPhenotypes();
};



class MessagesUtilities
{
public:
	static std::string generateStartingSentence(const std::vector<std::string>& tasks);
};
