#include "StatisticsGenerator.h"
#include <algorithm>



void StatisticsGenerator::fillTaskData()
{
	tasksData.reserve(params.stage1Params.GetNBins());
	for (uint32_t i = 0; i < params.stage1Params.GetNBins(); ++i)
	{
		tasksData.push_back(TaskData{ i });
	}
	std::vector<uint64_t> nOutputKmersPerBin;
	readDump(nOutputKmersPerBin, params.statisticsParams.statsNOutputKmers);
	std::sort(tasksData.begin(), tasksData.end(), [&](const TaskData& a, const TaskData& b) { return nOutputKmersPerBin[a.binId] > nOutputKmersPerBin[b.binId]; });
}



void StatisticsGenerator::operator()()
{
	TaskData taskData;
	while (tasksPool.getTask(taskData))
	{
		std::ifstream matrixFile(params.mkmcParams.outputMatrixFiles[taskData.binId]);
		if (!matrixFile.is_open())
		{
			std::cerr << "Error: cannot open " << params.mkmcParams.outputMatrixFiles[taskData.binId] << "." << std::endl;
			exit(1);
		}

		std::ofstream normFile(params.mkmcParams.outputFilesNorm[taskData.binId]);
		if (!normFile.is_open())
		{
			std::cerr << "Error: cannot open " << params.mkmcParams.outputFilesNorm[taskData.binId] << "." << std::endl;
			exit(1);
		}

		bool generatePearson = false, generateSpearman = false, generateKendall = false;
		bool generateEntropy = false;
		bool generateStatistics = false, generateTTest = false, generateSNR = false, generateWilcoxonRankSum = false, generateDIDS = false, generateANOVA = false;

		for (auto method : params.statisticsParams.correlationMethods)
		{
			if (method == StatisticsParams::CorrelationMethod::Pearson)
				generatePearson = true;
			else if (method == StatisticsParams::CorrelationMethod::Spearman)
				generateSpearman = true;
			else if (method == StatisticsParams::CorrelationMethod::Kendall)
				generateKendall = true;
		}
		generateEntropy = params.statisticsParams.generateEntropy;
		generateStatistics = !params.statisticsParams.classificationMethods.empty();
		for (auto method : params.statisticsParams.classificationMethods)
		{
			if (method == StatisticsParams::DifferentialAnalysisMethod::TTest)
				generateTTest = true;
			else if (method == StatisticsParams::DifferentialAnalysisMethod::SNR)
				generateSNR = true;
			else if (method == StatisticsParams::DifferentialAnalysisMethod::WilcoxonRankSum)
				generateWilcoxonRankSum = true;
			else if (method == StatisticsParams::DifferentialAnalysisMethod::DIDS)
				generateDIDS = true;
			else if (method == StatisticsParams::DifferentialAnalysisMethod::ANOVA)
				generateANOVA = true;
		}

		std::string header;
		std::getline(matrixFile, header);
		normFile << header << '\n';

		std::ofstream pearsonFile, spearmanFile, kendallFile;
		std::ofstream entropyFile;
		std::ofstream tTestFile, SNRFile, wilcoxonRankSumFile, DIDSFile, ANOVAFile;
		if (generatePearson)
		{
			pearsonFile.open(params.mkmcParams.outputFilesPearson[taskData.binId]);
			if (!pearsonFile.is_open())
			{
				std::cerr << "Error: cannot open " << params.mkmcParams.outputFilesPearson[taskData.binId] << "." << std::endl;
				exit(1);
			}
			pearsonFile << "k-mer\tcorrelation\n";
		}
		if (generateSpearman)
		{
			spearmanFile.open(params.mkmcParams.outputFilesSpearman[taskData.binId]);
			if (!spearmanFile.is_open())
			{
				std::cerr << "Error: cannot open " << params.mkmcParams.outputFilesSpearman[taskData.binId] << "." << std::endl;
				exit(1);
			}
			spearmanFile << "k-mer\tcorrelation\n";
		}
		if (generateKendall)
		{
			kendallFile.open(params.mkmcParams.outputFilesKendall[taskData.binId]);
			if (!kendallFile.is_open())
			{
				std::cerr << "Error: cannot open " << params.mkmcParams.outputFilesKendall[taskData.binId] << "." << std::endl;
				exit(1);
			}
			kendallFile << "k-mer\tcorrelation\n";
		}
		if (generateEntropy)
		{
			entropyFile.open(params.mkmcParams.outputFilesEntropy[taskData.binId]);
			if (!entropyFile.is_open())
			{
				std::cerr << "Error: cannot open " << params.mkmcParams.outputFilesEntropy[taskData.binId] << "." << std::endl;
				exit(1);
			}
			kendallFile << "k-mer\tentropy\n";
		}
		if (generateStatistics)
		{
			if (generateTTest)
			{
				tTestFile.open(params.mkmcParams.outputFilesTTest[taskData.binId]);
				if (!tTestFile.is_open())
				{
					std::cerr << "Error: cannot open " << params.mkmcParams.outputFilesTTest[taskData.binId] << "." << std::endl;
					exit(1);
				}
				tTestFile << "k-mer\tp-value\n";
			}
			if (generateSNR)
			{
				SNRFile.open(params.mkmcParams.outputFilesSNR[taskData.binId]);
				if (!SNRFile.is_open())
				{
					std::cerr << "Error: cannot open " << params.mkmcParams.outputFilesSNR[taskData.binId] << "." << std::endl;
					exit(1);
				}
				SNRFile << "k-mer\Signal to Noise ratio\n";
			}
			if (generateWilcoxonRankSum)
			{
				wilcoxonRankSumFile.open(params.mkmcParams.outputFilesWilcoxonRankSum[taskData.binId]);
				if (!wilcoxonRankSumFile.is_open())
				{
					std::cerr << "Error: cannot open " << params.mkmcParams.outputFilesWilcoxonRankSum[taskData.binId] << "." << std::endl;
					exit(1);
				}
				wilcoxonRankSumFile << "k-mer\tp-value\n";
			}
			if (generateDIDS)
			{
				DIDSFile.open(params.mkmcParams.outputFilesDIDS[taskData.binId]);
				if (!DIDSFile.is_open())
				{
					std::cerr << "Error: cannot open " << params.mkmcParams.outputFilesDIDS[taskData.binId] << "." << std::endl;
					exit(1);
				}
				DIDSFile << "k-mer\tDIDS\n";
			}
			if (generateANOVA)
			{
				ANOVAFile.open(params.mkmcParams.outputFilesANOVA[taskData.binId]);
				if (!ANOVAFile.is_open())
				{
					std::cerr << "Error: cannot open " << params.mkmcParams.outputFilesANOVA[taskData.binId] << "." << std::endl;
					exit(1);
				}
				ANOVAFile << "k-mer\tp-value\n";
			}
		}

		refresh::normalization_work<uint64_t, double> normalization;
		if (params.statisticsParams.generateNormalization)
		{
			normalization.register_method(params.statisticsParams.normalizationMethod);
			normalization.set_no_series(params.mkmcParams.samples.size());
			normalization.deserialize(params.statisticsParams.normalizationMethod, normalizationData);

			normalization.initialize();
		}

		refresh::correlation correlation;
		refresh::statistics_entropy entropyObj;
		refresh::statistical_test statistics;
		refresh::scorers scorer;

		std::string kmerSequence;
		std::vector<uint64_t> matrixEntry;
		std::vector<double> normEntry;
		matrixEntry.resize(params.mkmcParams.samples.size());
		normEntry.resize(params.mkmcParams.samples.size());

		ProgressBarUpdater progress_bar_updater(progress_bar, (std::max)(1ull, totAllKmers / 100ull));

		while (true)
		{
			if (!getLine(matrixFile, kmerSequence, matrixEntry))
				break;

			normalization.norm_entry(params.statisticsParams.normalizationMethod, matrixEntry, normEntry);
			putLine(normFile, kmerSequence, normEntry);

			if (generatePearson)
			{
				const double pearson = refresh::correlation::pearson(normEntry.begin(), normEntry.end(), correlationPhenotype.begin());
				putLine(pearsonFile, kmerSequence, { pearson });
			}
			if (generateSpearman)
			{
				const double spearman = correlation.spearman(normEntry.begin(), normEntry.end(), correlationPhenotype.begin());
				putLine(spearmanFile, kmerSequence, { spearman });
			}
			if (generateKendall)
			{
				const double kendall = refresh::correlation::kendall_tau(normEntry.begin(), normEntry.end(), correlationPhenotype.begin());
				putLine(kendallFile, kmerSequence, { kendall });
			}
			if (generateEntropy)
			{
				const double entropy = entropyObj.entropy(matrixEntry.begin(), matrixEntry.end());
				putLine(entropyFile, kmerSequence, { entropy });
			}
			if (generateStatistics)
			{
				if (generateTTest)
				{
					const double tTestPValue = statistics.t_test(matrixEntry.begin(), matrixEntry.end(), differentialAnalysisPhenotype.begin()).p_value;
					putLine(tTestFile, kmerSequence, { tTestPValue });
				}
				if (generateSNR)
				{
					const double SNRPValue = statistics.SNR_test(matrixEntry.begin(), matrixEntry.end(), differentialAnalysisPhenotype.begin());
					putLine(SNRFile, kmerSequence, { SNRPValue });
				}
				if (generateWilcoxonRankSum)
				{
					const double wilcoxonRankSumPValue = statistics.mann_whitney_U_test(matrixEntry.begin(), matrixEntry.end(), differentialAnalysisPhenotype.begin()).p_value;
					putLine(wilcoxonRankSumFile, kmerSequence, { wilcoxonRankSumPValue });
				}
				if (generateDIDS)
				{
					const double dids = scorer.dids(matrixEntry.begin(), matrixEntry.end(), differentialAnalysisPhenotype.begin(), differentialAnalysisClasses);
					putLine(DIDSFile, kmerSequence, { dids });
				}
				if (generateANOVA)
				{
					const double anova = scorer.anova(matrixEntry.begin(), matrixEntry.end(), differentialAnalysisPhenotype.begin(), differentialAnalysisClasses).p_value;
					putLine(ANOVAFile, kmerSequence, { anova });
				}
			}

			++progress_bar_updater;
		}
	}
}


void StatisticsGenerator::generateStatisticsParallel()
{
	fillTaskData();

	if (params.statisticsParams.generateNormalization)
	{
		if (params.statisticsParams.normalizationMethod == StatisticsParams::NormalizationMethod::frequency_count)
			readDump(normalizationData, params.statisticsParams.normFrequencyFileTmp);
		else if (params.statisticsParams.normalizationMethod == StatisticsParams::NormalizationMethod::quantile)
			readDump(normalizationData, params.statisticsParams.normQuantileFileTmp);
	}

	std::vector<std::thread> threads(params.mkmcParams.nThreads);
	for (uint32_t i_thred = 0; i_thred < params.mkmcParams.nThreads; ++i_thred)
	{
		threads[i_thred] = std::thread([this] { (*this)(); });
	}

	for (std::thread& thread : threads)
	{
		thread.join();
	}
}
