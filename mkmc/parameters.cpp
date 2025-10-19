#include "parameters.h"
#include "SamplesFileReader.h"
#include "refresh/deterministic_random/lib/deterministic_random.h"
#include "Logger.h"
#include <iostream>
#include <filesystem>
#include <algorithm>
#include <thread>



MKMCParams::MKMCParams() :
	nThreads((std::min)(16U, std::thread::hardware_concurrency()))
{}



void StatisticsParams::CVParams::generateSamplesToExcludeOrder(const size_t nSamples)
{
	samplesToExcludeOrder.resize(nSamples);
	std::generate(samplesToExcludeOrder.begin(), samplesToExcludeOrder.end(), []() { static size_t n = 0; return n++; });
	if (p != 1) // For LOOCV random order is not necessary
	{
		std::mt19937 gen(seed);
		partial_shuffle(samplesToExcludeOrder.begin(), samplesToExcludeOrder.end(), samplesToExcludeOrder.end(), gen); // deterministic portable random
	}
}



std::string StatisticsParams::CVParams::getOutputFileNameImpl(CorrelationMethod method, size_t nSamples, size_t iTest, size_t nTests) const
{
	std::string methodStr;
	switch (method)
	{
	case CorrelationMethod::Pearson: methodStr = "pearson"; break;
	case CorrelationMethod::Spearman: methodStr = "spearman"; break;
	case CorrelationMethod::Kendall: methodStr = "kendall_tau"; break;
	}

	return outputFilesTemplate + "_cv_" + methodStr + "_" + std::to_string(iTest + 1) + "_" + std::to_string(nTests);
}



Params::Params() :
	phenotypes({*this,*this})
{
	stage1Params.SetReopenTmeEachTime(true);

	stage1Params.SetSignatureSelectionScheme(KMC::SignatureSelectionScheme::min_hash);

	// set KMC defaults
	stage1Params.SetCutoffMin(defaultKMCParams.ci);
	stage1Params.SetCutoffMax(defaultKMCParams.cx);
	stage1Params.SetCounterMax(defaultKMCParams.cs);
	stage1Params.SetKmerLen(defaultKMCParams.k);

	static KMC::NullPercentProgressObserver nullPercentProgressObserver;
	static KMC::NullProgressObserver nullProgressObserver;
	stage1Params.SetPercentProgressObserver(&nullPercentProgressObserver);
	stage1Params.SetProgressObserver(&nullProgressObserver);
}



bool Params::readAdditionalDataFromFiles(bool& warningPrinted)
{
	SamplesFileReader tasksFiller(mkmcParams);
	return tasksFiller.readSamples(mkmcParams.samples, warningPrinted);
}



void Params::generateTempAndOutputFilesNames()
{
	std::string tmpFilesTemplate = mkmcParams.tmpPath;
	if (mkmcParams.tmpPath.back() != '/' && mkmcParams.tmpPath.back() != '\\')
	{
		tmpFilesTemplate += static_cast<char>(std::filesystem::path::preferred_separator);
	}

	for (uint32_t tmp_database_id = 0; tmp_database_id < mkmcParams.samples.size(); ++tmp_database_id)
	{
		std::ostringstream sstreamKMCDir, sstreamKMC;
		sstreamKMCDir << tmpFilesTemplate << "kmc_tmp_" << std::setfill('0') << std::setw(5) << tmp_database_id;
		sstreamKMC << tmpFilesTemplate << "kmc_db_" << std::setfill('0') << std::setw(5) << tmp_database_id;

		mkmcParams.kmcTmpDirs.push_back(sstreamKMCDir.str());
		mkmcParams.kmcOutputFiles.push_back(sstreamKMC.str());
	}

	filterParams.kmersSequencesToFilterOutDB = mkmcParams.tmpPath + static_cast<char>(std::filesystem::path::preferred_separator) + "filter";
	mutableParams.kmersSequencesToFilterOut = mkmcParams.tmpPath + static_cast<char>(std::filesystem::path::preferred_separator) + "filter.fa";

	mkmcParams.outputMatrixBinFile = mkmcParams.outputFilesTemplate + ".kmcdb";
	mkmcParams.outputMatrixFile = mkmcParams.outputFilesTemplate + "_matrix";
	mkmcParams.outputFASTAFile = mkmcParams.outputFilesTemplate + ".fa";

	mkmcParams.normLearningBinFile = mkmcParams.outputFilesTemplate + ".stats";
	mkmcParams.normLearningBinFileSupplemented = mkmcParams.outputFilesTemplate + ".stats_supp";

	mkmcParams.outputFileNorm = mkmcParams.outputFilesTemplate + "_matrix_norm";

	mkmcParams.outputFilePearson = mkmcParams.outputFilesTemplate + "_pearson";
	mkmcParams.outputFilePearsonTop = mkmcParams.outputFilesTemplate + "_pearson_top";
	mkmcParams.outputFilePearsonTopCntMatrix = mkmcParams.outputFilesTemplate + "_pearson_top_matrix";
	mkmcParams.outputFilePearsonTopFasta = mkmcParams.outputFilesTemplate + "_pearson_top.fa";

	mkmcParams.outputFileSpearman = mkmcParams.outputFilesTemplate + "_spearman";
	mkmcParams.outputFileSpearmanTop = mkmcParams.outputFilesTemplate + "_spearman_top";
	mkmcParams.outputFileSpearmanTopCntMatrix = mkmcParams.outputFilesTemplate + "_spearman_top_matrix";
	mkmcParams.outputFileSpearmanTopFasta = mkmcParams.outputFilesTemplate + "_spearman_top.fa";

	mkmcParams.outputFileKendall = mkmcParams.outputFilesTemplate + "_kendall_tau";
	mkmcParams.outputFileKendallTop = mkmcParams.outputFilesTemplate + "_kendall_tau_top";
	mkmcParams.outputFileKendallTopCntMatrix = mkmcParams.outputFilesTemplate + "_kendall_tau_top_matrix";
	mkmcParams.outputFileKendallTopFasta = mkmcParams.outputFilesTemplate + "_kendall_tau_top.fa";


	mkmcParams.outputFileEntropy = mkmcParams.outputFilesTemplate + "_entropy";
	mkmcParams.outputFileEntropyTop = mkmcParams.outputFilesTemplate + "_entropy_top";
	mkmcParams.outputFileEntropyTopCntMatrix = mkmcParams.outputFilesTemplate + "_entropy_top_matrix";
	mkmcParams.outputFileEntropyTopFasta = mkmcParams.outputFilesTemplate + "_entropy_top.fa";

	mkmcParams.outputFileTTest = mkmcParams.outputFilesTemplate + "_ttest";
	mkmcParams.outputFileTTestCor = mkmcParams.outputFilesTemplate + "_ttest_cor_all";
	mkmcParams.outputFileTTestCorSignificant = mkmcParams.outputFilesTemplate + "_ttest_cor_significant";
	mkmcParams.outputFileTTestCorSignificantCntMatrix = mkmcParams.outputFilesTemplate + "_ttest_cor_significant_matrix";
	mkmcParams.outputFileTTestCorSignificantFasta = mkmcParams.outputFilesTemplate + "_ttest_cor_significant.fa";

	mkmcParams.outputFileSNR = mkmcParams.outputFilesTemplate + "_snr";
	mkmcParams.outputFileSNRTop = mkmcParams.outputFilesTemplate + "_snr_top";
	mkmcParams.outputFileSNRTopCntMatrix = mkmcParams.outputFilesTemplate + "_snr_top_matrix";
	mkmcParams.outputFileSNRTopFasta = mkmcParams.outputFilesTemplate + "_snr_top.fa";

	mkmcParams.outputFileUnnormalizedSNR = mkmcParams.outputFilesTemplate + "_snr_for_unnornalized";
	mkmcParams.outputFileUnnormalizedSNRTop = mkmcParams.outputFilesTemplate + "_snr_top_for_unnornalized";
	mkmcParams.outputFileUnnormalizedSNRTopCntMatrix = mkmcParams.outputFilesTemplate + "_snr_top_matrix_for_unnornalized";
	mkmcParams.outputFileUnnormalizedSNRTopFasta = mkmcParams.outputFilesTemplate + "_snr_top_for_unnornalized.fa";

	mkmcParams.outputFileWilcoxonRankSum = mkmcParams.outputFilesTemplate + "_wrs";
	mkmcParams.outputFileWilcoxonRankSumCor = mkmcParams.outputFilesTemplate + "_wrs_cor_all";
	mkmcParams.outputFileWilcoxonRankSumCorSignificant = mkmcParams.outputFilesTemplate + "_wrs_cor_significant";
	mkmcParams.outputFileWilcoxonRankSumCorSignificantCntMatrix = mkmcParams.outputFilesTemplate + "_wrs_cor_significant_matrix";
	mkmcParams.outputFileWilcoxonRankSumCorSignificantFasta = mkmcParams.outputFilesTemplate + "_wrs_cor_significant.fa";

	mkmcParams.outputFileDIDS = mkmcParams.outputFilesTemplate + "_dids";
	mkmcParams.outputFileDIDSTop = mkmcParams.outputFilesTemplate + "_dids_top";
	mkmcParams.outputFileDIDSTopCntMatrix = mkmcParams.outputFilesTemplate + "_dids_top_matrix";
	mkmcParams.outputFileDIDSTopFasta = mkmcParams.outputFilesTemplate + "_dids_top.fa";

	mkmcParams.outputFileANOVA = mkmcParams.outputFilesTemplate + "_anova";
	mkmcParams.outputFileANOVACor = mkmcParams.outputFilesTemplate + "_anova_cor_all";
	mkmcParams.outputFileANOVACorSignificant = mkmcParams.outputFilesTemplate + "_anova_cor_significant";
	mkmcParams.outputFileANOVACorSignificantCntMatrix = mkmcParams.outputFilesTemplate + "_anova_cor_significant_matrix";
	mkmcParams.outputFileANOVACorSignificantFasta = mkmcParams.outputFilesTemplate + "_anova_cor_significant.fa";

	mkmcParams.outputFileUMAP = mkmcParams.outputFilesTemplate + "_umap";
	mkmcParams.outputFilePCA = mkmcParams.outputFilesTemplate + "_pca";
	mkmcParams.outputFilePCAVariance = mkmcParams.outputFilesTemplate + "_pca_variance";

	mkmcParams.outputFileTotCnt = mkmcParams.outputFilesTemplate + "_tot_cnt";
}



bool Params::adjustKMCPerformanceParams()
{
	bool warningPrinted = false;
	bool mKMCWorkersReduced = false;
	if (mkmcParams.nThreads == 1)
	{
		mkmcParams.nKMCWorkers = 1;
		mKMCWorkersReduced = true;
	}
	else if (mkmcParams.nKMCWorkers > mkmcParams.nThreads - 1)
	{
		mkmcParams.nKMCWorkers = mkmcParams.nThreads - 1;
		mKMCWorkersReduced = true;
	}

	stage1Params.SetNThreads(mkmcParams.nThreads / mkmcParams.nKMCWorkers);
	stage2Params.SetNThreads(mkmcParams.nThreads / mkmcParams.nKMCWorkers);

	if (mkmcParams.nKMCWorkers * 2 > mkmcParams.maxRamGB)
	{
		mkmcParams.nKMCWorkers = mkmcParams.maxRamGB / 2;
		mKMCWorkersReduced = true;
	}

	if (mkmcParams.nKMCWorkersUserSet && mKMCWorkersReduced)
	{
		Logger::Inst().Log("Warning: number of workers is too huge, reduced to " + std::to_string(mkmcParams.nKMCWorkers) + ".", 1);
		warningPrinted = true;
	}

	stage1Params.SetMaxRamGB(mkmcParams.maxRamGB / mkmcParams.nKMCWorkers);
	stage2Params.SetMaxRamGB(mkmcParams.maxRamGB / mkmcParams.nKMCWorkers);

	stage1Params.SetNBins(mkmcParams.nKMCBins);

	return warningPrinted;
}



bool Params::adjustAnotherParams()
{
	bool warningPrinted = false;
	if (mkmcParams.maxRamGBUserDefined && stage1Params.GetRamOnlyMode())
	{
		Logger::Inst().Log("Warning: when -r parameter is given, the limit specified with -m may be exceeded.", 1);
		warningPrinted = true;
	}

	size_t nCorrelationMethods = statisticsParams.correlationMethods.size();
	std::sort(statisticsParams.correlationMethods.begin(), statisticsParams.correlationMethods.end());
	statisticsParams.correlationMethods.erase(std::unique(statisticsParams.correlationMethods.begin(), statisticsParams.correlationMethods.end()), statisticsParams.correlationMethods.end());
	if (nCorrelationMethods != statisticsParams.correlationMethods.size())
	{
		Logger::Inst().Log("Warning: some correlation methods were given multiple times.", 1);
		warningPrinted = true;
	}

	size_t nOutputFileTypes = mkmcParams.outputFileTypes.size();
	std::sort(mkmcParams.outputFileTypes.begin(), mkmcParams.outputFileTypes.end());
	mkmcParams.outputFileTypes.erase(std::unique(mkmcParams.outputFileTypes.begin(), mkmcParams.outputFileTypes.end()), mkmcParams.outputFileTypes.end());

	if (nOutputFileTypes != mkmcParams.outputFileTypes.size())
	{
		Logger::Inst().Log("Warning: some output files types were given multiple times.", 1);
		warningPrinted = true;
	}

	statisticsParams.umap_params.num_threads = mkmcParams.nThreads;

	statisticsParams.cvParams.generateSamplesToExcludeOrder(mkmcParams.samples.size());

	return warningPrinted;
}



bool Params::readPhenotypes()
{
	if (!phenotypes.correlationPhenotype.getFileName().empty())
		if (!phenotypes.correlationPhenotype.readPhenotype())
			return false;

	if (!phenotypes.differentialAnalysisPhenotype.getFileName().empty())
	{
		if (!phenotypes.differentialAnalysisPhenotype.readPhenotype())
			return false;
		if (!phenotypes.differentialAnalysisPhenotype.mapPhenotypeToInts())
			return false;

		for (auto method : statisticsParams.classificationMethods)
		{
			if (method == StatisticsParams::DifferentialAnalysisMethod::TTest && phenotypes.differentialAnalysisPhenotype.getClassesNumber() > 2)
			{
				Logger::Inst().Log("Error: number of distinct classes in a file " + phenotypes.differentialAnalysisPhenotype.getFileName() + " for T-Test must equal to 2.");
				return false;
			}
			else if (method == StatisticsParams::DifferentialAnalysisMethod::SNR && phenotypes.differentialAnalysisPhenotype.getClassesNumber() > 2)
			{
				Logger::Inst().Log("Error: number of distinct classes in a file " + phenotypes.differentialAnalysisPhenotype.getFileName() + " for Signal to Noise ratio determination must equal to 2.");
				return false;
			}
			else if (method == StatisticsParams::DifferentialAnalysisMethod::WilcoxonRankSum && phenotypes.differentialAnalysisPhenotype.getClassesNumber() > 2)
			{
				Logger::Inst().Log("Error: number of distinct classes in a file " + phenotypes.differentialAnalysisPhenotype.getFileName() + " for Wilcoxon-rank sum determination must equal to 2.");
				return false;
			}
		}
	}

	return true;
}



std::string MessagesUtilities::generateSentence(const std::vector<std::string>& tasks, bool capitalize/* = false*/)
{
	std::string result;
	if (tasks.size() == 1)
		result += *tasks.begin();
	else if (tasks.size() == 2)
		result += *tasks.begin() + " and " + result += *(tasks.begin() + 1);
	else
	{
		for (size_t i = 0; i < tasks.size(); ++i)
		{
			if (i == tasks.size() - 1)
				result += ", and ";
			else if (i != 0)
				result += ", ";
			result += tasks[i];
		}
	}
	if (capitalize)
		result.front() = std::toupper(result.front());
	return result;
}

