#include "parameters.h"
#include "SamplesFileReader.h"
#include <iostream>
#include <filesystem>
#include <algorithm>



Params::Params() :
	phenotypes({*this,*this})
{
	stage1Params.SetReopenTmeEachTime(true);

	stage1Params.SetSignatureSelectionScheme(KMC::SignatureSelectionScheme::min_hash);

	stage1Params.SetCutoffMin(1);
	stage1Params.SetCutoffMax(static_cast<uint64_t>(4E9));
	stage1Params.SetCounterMax(65535);

	static KMC::NullPercentProgressObserver nullPercentProgressObserver;
	static KMC::NullProgressObserver nullProgressObserver;
	stage1Params.SetPercentProgressObserver(&nullPercentProgressObserver);
	stage1Params.SetProgressObserver(&nullProgressObserver);
}



bool Params::readAdditionalParamsFromFiles()
{
	SamplesFileReader tasksFiller(mkmcParams);
	return tasksFiller.readSamples(mkmcParams.samples);
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

	mkmcParams.outputBinFile = mkmcParams.outputFilesTemplate + ".kmcdb";
	mkmcParams.outputStatsBinFile = mkmcParams.outputFilesTemplate + "_norm+cor.kmcdb"; 
	mkmcParams.outputMatrixFile = mkmcParams.outputFilesTemplate + "_matrix";
	mkmcParams.outputFASTAFile = mkmcParams.outputFilesTemplate + ".fa";

	mkmcParams.outputFileNorm = mkmcParams.outputFilesTemplate + "_norm";

	mkmcParams.outputFilePearson = mkmcParams.outputFilesTemplate + "_pearson";
	mkmcParams.outputFileSpearman = mkmcParams.outputFilesTemplate + "_spearman";
	mkmcParams.outputFileKendall = mkmcParams.outputFilesTemplate + "_kendall_tau";

	mkmcParams.outputFileEntropy = mkmcParams.outputFilesTemplate + "_entropy";

	mkmcParams.outputFileTTest = mkmcParams.outputFilesTemplate + "_ttest";
	mkmcParams.outputFileSNR = mkmcParams.outputFilesTemplate + "_snr";
	mkmcParams.outputFileWilcoxonRankSum = mkmcParams.outputFilesTemplate + "_wrs";
	mkmcParams.outputFileDIDS = mkmcParams.outputFilesTemplate + "_dids";
	mkmcParams.outputFileANOVA = mkmcParams.outputFilesTemplate + "_anova";
}



void Params::adjustKMCPerformanceParams()
{
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
		std::cerr << "Warning: number of workers is too huge, reduced to " << mkmcParams.nKMCWorkers << "." << std::endl;
	}

	stage1Params.SetMaxRamGB(mkmcParams.maxRamGB / mkmcParams.nKMCWorkers);
	stage2Params.SetMaxRamGB(mkmcParams.maxRamGB / mkmcParams.nKMCWorkers);

	stage1Params.SetNBins(mkmcParams.nKMCBins);
}



void Params::adjustAnotherParams()
{
	if (mkmcParams.maxRamGBUserDefined && stage1Params.GetRamOnlyMode())
	{
		std::cerr << "Warning: when -r parameter is given, limit specified with -m may be exceeded." << std::endl;
	}

	size_t nCorrelationMethods = statisticsParams.correlationMethods.size();
	std::sort(statisticsParams.correlationMethods.begin(), statisticsParams.correlationMethods.end());
	statisticsParams.correlationMethods.erase(std::unique(statisticsParams.correlationMethods.begin(), statisticsParams.correlationMethods.end()), statisticsParams.correlationMethods.end());
	if (nCorrelationMethods != statisticsParams.correlationMethods.size())
	{
		std::cerr << "Warning: some correlation methods were given multiple times." << std::endl;
	}

	size_t nOutputFileTypes = mkmcParams.outputFileTypes.size();
	std::sort(mkmcParams.outputFileTypes.begin(), mkmcParams.outputFileTypes.end());
	mkmcParams.outputFileTypes.erase(std::unique(mkmcParams.outputFileTypes.begin(), mkmcParams.outputFileTypes.end()), mkmcParams.outputFileTypes.end());

	if (nOutputFileTypes != mkmcParams.outputFileTypes.size())
	{
		std::cerr << "Warning: some output files types were given multiple times." << std::endl;
	}

	if (statisticsParams.generateNormalization)
	{
		if (std::find(mkmcParams.outputFileTypes.begin(), mkmcParams.outputFileTypes.end(), OutputFileType::Matrix) == mkmcParams.outputFileTypes.end())
		{
			std::cerr << "Warning: due to normalization generation, temporarily MKMC has to generate output matrix (-o matrix flag will be additionally applied)." << std::endl;
			mkmcParams.outputFileTypes.push_back(OutputFileType::Matrix);
		}
	}
}



void Params::readPhenotypes()
{
	if (!phenotypes.correlationPhenotype.getFileName().empty())
		phenotypes.correlationPhenotype.readPhenotype();

	if (!phenotypes.differentialAnalysisPhenotype.getFileName().empty())
	{
		phenotypes.differentialAnalysisPhenotype.readPhenotype();
		phenotypes.differentialAnalysisPhenotype.mapPhenotypeToInts();

		for (auto method : statisticsParams.classificationMethods)
		{
			if (method == StatisticsParams::DifferentialAnalysisMethod::TTest && phenotypes.differentialAnalysisPhenotype.getClassesNumber() > 2)
			{
				std::cerr << "Error: number of distinct classes in a file " << phenotypes.differentialAnalysisPhenotype.getFileName() << " for T-Test must equal to 2." << std::endl;
				exit(1);
			}
			else if (method == StatisticsParams::DifferentialAnalysisMethod::SNR && phenotypes.differentialAnalysisPhenotype.getClassesNumber() > 2)
			{
				std::cerr << "Error: number of distinct classes in a file " << phenotypes.differentialAnalysisPhenotype.getFileName() << " for Signal to Noise ratio determination must equal to 2." << std::endl;
				exit(1);
			}
			else if (method == StatisticsParams::DifferentialAnalysisMethod::WilcoxonRankSum && phenotypes.differentialAnalysisPhenotype.getClassesNumber() > 2)
			{
				std::cerr << "Error: number of distinct classes in a file " << phenotypes.differentialAnalysisPhenotype.getFileName() << " for Wilcoxon-rank sum determination must equal to 2." << std::endl;
				exit(1);
			}
		}
	}
}
