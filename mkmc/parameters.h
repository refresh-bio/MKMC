#pragma once
#include <vector>
#include <string>
#include <cstdint>
#include "PhenotypeReaders.h"
#include "kmc_core/kmc_runner.h"
#include "kmc_api/kmc_file.h"
#include "kmc_api/kmer_api.h"
#undef small
#include "refresh/statistics/lib/statistics_normalization.h"
#include "refresh/statistics/lib/statistics_umap.h"
#include "refresh/statistics/lib/statistics_pca.h"


enum class OutputFileType { FASTA, Matrix, _N }; // _N - number of possibilities

struct Sample
{
	std::string name;
	std::vector<std::string> inputFiles;
};

struct MKMCParams
{
	std::string inputFileName;
	std::string outputFilesTemplate;
	std::string tmpPath;

	std::string getLogFileName() const {
		return outputFilesTemplate + ".log";
	}

	std::vector<Sample> samples;

	KMC::InputFileType inputFileType = KMC::InputFileType::FASTQ; // also default "fq" value in input parameters
	std::vector<std::string> kmcOutputFiles;
	std::vector<std::string> kmcTmpDirs;

	std::string outputMatrixBinFile;
	std::string outputFASTAFile;
	std::string outputMatrixFile;
	std::vector<OutputFileType> outputFileTypes;

	std::string normLearningBinFile;
	std::string normLearningBinFileSupplemented;

	std::string outputFileNorm;

	std::string outputFilePearson;
	std::string outputFilePearsonTop;
	std::string outputFilePearsonTopCntMatrix;
	std::string outputFilePearsonTopFasta;

	std::string outputFileSpearman;
	std::string outputFileSpearmanTop;
	std::string outputFileSpearmanTopCntMatrix;
	std::string outputFileSpearmanTopFasta;

	std::string outputFileKendall;
	std::string outputFileKendallTop;
	std::string outputFileKendallTopCntMatrix;
	std::string outputFileKendallTopFasta;

	std::string outputFileEntropy;
	std::string outputFileEntropyTop;
	std::string outputFileEntropyTopCntMatrix;
	std::string outputFileEntropyTopFasta;

	std::string outputFileTTest;
	std::string outputFileTTestCor;
	std::string outputFileTTestCorSignificant;
	std::string outputFileTTestCorSignificantCntMatrix;
	std::string outputFileTTestCorSignificantFasta;

	std::string outputFileSNR;
	std::string outputFileSNRTop;
	std::string outputFileSNRTopCntMatrix;
	std::string outputFileSNRTopFasta;

	std::string outputFileUnnormalizedSNR;
	std::string outputFileUnnormalizedSNRTop;
	std::string outputFileUnnormalizedSNRTopCntMatrix;
	std::string outputFileUnnormalizedSNRTopFasta;

	std::string outputFileWilcoxonRankSum;
	std::string outputFileWilcoxonRankSumCor;
	std::string outputFileWilcoxonRankSumCorSignificant;
	std::string outputFileWilcoxonRankSumCorSignificantCntMatrix;
	std::string outputFileWilcoxonRankSumCorSignificantFasta;

	std::string outputFileDIDS;
	std::string outputFileDIDSTop;
	std::string outputFileDIDSTopCntMatrix;
	std::string outputFileDIDSTopFasta;

	std::string outputFileANOVA;
	std::string outputFileANOVACor;
	std::string outputFileANOVACorSignificant;
	std::string outputFileANOVACorSignificantCntMatrix;
	std::string outputFileANOVACorSignificantFasta;

	std::string outputFileUMAP;
	std::string outputFilePCA;
	std::string outputFilePCAVariance;

	bool totCntGeneration = false;
	std::string outputFileTotCnt;

	uint32_t nThreads; // initialized in constructor
	uint32_t nKMCWorkers = 4;
	uint32_t maxRamGB = 16;
	bool maxRamGBUserDefined = false;

	uint32_t nKMCBins = 512;

	bool nKMCWorkersUserSet = false;

	bool reuseDBFiles = false;
	bool keepKMCdbs = false;

	bool generateForNonNormalized = false;

	int verbosity_level = 0;

	MKMCParams();
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

	static std::string getNormalizationMethodStreamName(const NormalizationMethod normalizationMethod)
	{
		switch (normalizationMethod) {
		case NormalizationMethod::frequency_count: return "norm_frequency";
		case NormalizationMethod::quantile: return "norm_quantile";
		case NormalizationMethod::deseq2: return "norm_deseq2";
		default: assert(false);
		}
		return "";
	}

	static std::vector<NormalizationMethod> getAllSupportedNormalizationMethods()
	{
		return { NormalizationMethod::frequency_count, NormalizationMethod::quantile, NormalizationMethod::deseq2 };
	}

	static std::vector<NormalizationMethod> getAlwaysLearnedNormalizationMethods()
	{
		return { NormalizationMethod::frequency_count };
	}

	static std::string getNormalizationMethodMKMCParamName(const NormalizationMethod normalizationMethod)
	{
		switch (normalizationMethod) {
		case NormalizationMethod::frequency_count: return "freq";
		case NormalizationMethod::quantile: return "q";
		case NormalizationMethod::deseq2: return "deseq";
		default: assert(false);
		}
		return "";
	}

	bool generateNormalization = false;
	NormalizationMethod normalizationMethod = NormalizationMethod::frequency_count; // initialization due to compiler warnings
	bool saveNormalization = false;
	bool learnDeseq2 = false;
	bool learnQuantile = false;
	bool normalizationLearningWasSupplemented = false;

	bool runUMAP = false;
	bool runPCA = false;

	bool nDimensionReductionUserDefined = false;
	uint32_t nDimensionReduction = 2;

	refresh::umap<double>::params_t umap_params{};  // default initialize="spectral" value in input parameters
	refresh::pca<double>::computation_mode_t pca_mod = refresh::pca<double>::computation_mode_t::svd; // also default "svd" value in input parameters

	enum class CorrelationMethod { Pearson, Spearman, Kendall };
	std::vector<CorrelationMethod> correlationMethods;

	size_t nTop = 10000;
	bool nTopUserDefined = false;

	struct CVParams
	{
		bool cv = false;
		size_t nTestSamples = 2;

		uint64_t seed = 1234567890;
		bool seedUserDefined = false;

		std::vector<size_t> samplesToBeTestOrder;

		CVParams(const std::string& outputFilesTemplate) : outputFilesTemplate(outputFilesTemplate) {}

		std::string getOuputFileNameTop(CorrelationMethod method, size_t nSamples, size_t iTest, size_t nFolds) const
		{
			return getOutputFileNameImpl(method, nSamples, iTest, nFolds) + "_top";
		}
		std::string getOuputFileNameTopCntMatrix(CorrelationMethod method, size_t nSamples, size_t iTest, size_t nFolds) const
		{
			return getOutputFileNameImpl(method, nSamples, iTest, nFolds) + "_top_matrix";
		}
		std::string getOuputFileNameTopFasta(CorrelationMethod method, size_t nSamples, size_t iTest, size_t nFolds) const
		{
			return getOutputFileNameImpl(method, nSamples, iTest, nFolds) + "_top.fa";
		}

		void generateSamplesToExcludeOrder(const size_t nSamples);

	private:
		const std::string& outputFilesTemplate;

		std::string getOutputFileNameImpl(CorrelationMethod method, size_t nSamples, size_t iTest, size_t nFolds) const;
	} cvParams;

	enum class DifferentialAnalysisMethod { TTest, SNR, WilcoxonRankSum, DIDS, ANOVA };
	std::vector<DifferentialAnalysisMethod> classificationMethods;
	bool didsModeUserDefined = false;
	enum class DIDSMode { sqrt, quadratic, tanh };
	DIDSMode didsMode = DIDSMode::sqrt; // also default "sqrt" value in input parameters

	bool correctPvalues = false;
	enum class DifferentialAnalysisCorrectionMethod { Bonferroni, HolmBonferroni, BenjaminiHochberg, BenjaminiYekutieli };
	DifferentialAnalysisCorrectionMethod classificationPValueCorrection = DifferentialAnalysisCorrectionMethod::Bonferroni; // initialization due to compiler warnings
	double maxCorrectedPval = 0.05;

	bool generateEntropy = false;

	StatisticsParams(const std::string& outputFilesTemplate) : cvParams(outputFilesTemplate) {}
};

struct MutableParams
{
	bool tmpDirCreated = false;

	bool kmcDbsCreated = false;

	mutable bool createdFastaFile = false;
	mutable std::string kmersSequencesToFilterOut; // file containing k-mers to count (input or generated)
};

struct DefaultKMCParams
{
	const uint32_t k = 25;
	const uint64_t ci = 1;
	const uint64_t cx = static_cast<uint64_t>(4e9);
	const uint64_t cs = 65535;
};

struct Phenotypes
{
	PhenotypeReader<double> correlationPhenotype;
	DifferentialAnalysisPhenotypeReader differentialAnalysisPhenotype; 
};

struct Params
{
	MKMCParams mkmcParams;
	mutable MutableParams mutableParams;

	static inline DefaultKMCParams defaultKMCParams;

	KMC::Stage1Params stage1Params;
	KMC::Stage2Params stage2Params;
	uint32_t ci = static_cast<uint32_t>(defaultKMCParams.ci); // cx and cs are set by CLI11 in KMC params; as filtering is sensitive for ci, we treat it differently

	FilterParams filterParams;
	StatisticsParams statisticsParams = mkmcParams.outputFilesTemplate;

	Phenotypes phenotypes;

	Params();
	void generateTempAndOutputFilesNames();
	bool readAdditionalDataFromFiles(bool& warningPrinted);
	bool adjustKMCPerformanceParams();
	bool adjustAnotherParams();
	bool readPhenotypes();
};



class MessagesUtilities
{
public:
	static bool generateSentence(const std::vector<std::string>& tasks, std::string& result, bool capitalize = false);
};
