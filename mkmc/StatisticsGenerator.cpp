#include "StatisticsGenerator.h"



void StatisticsGenerator::readPhenotype(std::vector<int>& phenotype)
{
	std::ifstream phenotypeFile(params.statisticsParams.phenotypeFile);
	if (!phenotypeFile.is_open())
	{
		std::cerr << "Error: cannot open " << params.statisticsParams.phenotypeFile << "." << std::endl;
		exit(1);
	}

	int value;
	while (phenotypeFile >> value)
		phenotype.push_back(value);

	if (phenotype.size() != params.mkmcParams.inputFilesPerSample.size())
	{
		std::cerr << "Error: a phenotype size in a file " << params.statisticsParams.phenotypeFile  <<  " (" << phenotype.size() << ") is different than number of samples (" << params.mkmcParams.inputFilesPerSample.size() << ")." << std::endl;
		exit(1);
	}
}



void StatisticsGenerator::readNormalizationDump(std::vector<uint8_t>& normalizationData, std::string normalizationFileName)
{
	std::ifstream normalizationFile(normalizationFileName, std::ios::binary);
	if (!normalizationFile.is_open())
	{
		std::cerr << "Error: cannot open " << normalizationFileName << "." << std::endl;
		exit(1);
	}
	size_t normalizationDataSize;
	normalizationFile.read(reinterpret_cast<char*>(&normalizationDataSize), sizeof(size_t));
	normalizationData.resize(normalizationDataSize);
	normalizationFile.read(reinterpret_cast<char*>(normalizationData.data()), normalizationDataSize * sizeof(uint8_t));
}



void StatisticsGenerator::operator()()
{
	TaskData taskData;
	while (tasksPool.getTask(taskData))
	{
		std::ifstream matrixFile(params.mkmcParams.outputFiles[taskData.binId]);
		if (!matrixFile.is_open())
		{
			std::cerr << "Error: cannot open " << params.mkmcParams.outputFiles[taskData.binId] << "." << std::endl;
			exit(1);
		}

		std::ofstream frequenceNormFile(params.mkmcParams.outputFilesNormFrequency[taskData.binId]);
		if (!frequenceNormFile.is_open())
		{
			std::cerr << "Error: cannot open " << params.mkmcParams.outputFilesNormFrequency[taskData.binId] << "." << std::endl;
			exit(1);
		}
		std::ofstream pearsonFile(params.mkmcParams.outputFilesPearson[taskData.binId]);
		if (!pearsonFile.is_open())
		{
			std::cerr << "Error: cannot open " << params.mkmcParams.outputFilesPearson[taskData.binId] << "." << std::endl;
			exit(1);
		}
		std::ofstream spearmanFile(params.mkmcParams.outputFilesSpearman[taskData.binId]);
		if (!spearmanFile.is_open())
		{
			std::cerr << "Error: cannot open " << params.mkmcParams.outputFilesSpearman[taskData.binId] << "." << std::endl;
			exit(1);
		}
		std::ofstream kendallFile(params.mkmcParams.outputFilesKendall[taskData.binId]);
		if (!kendallFile.is_open())
		{
			std::cerr << "Error: cannot open " << params.mkmcParams.outputFilesKendall[taskData.binId] << "." << std::endl;
			exit(1);
		}

		std::string header;
		std::getline(matrixFile, header);
		frequenceNormFile << header << '\n';
		pearsonFile << header << '\n';

		refresh::normalization_work<size_t, double> normalization;

		normalization.register_method(StatisticsParams::NormalizationMethod::frequency_count);
		//norm.register_method(StatisticsParams::NormalizationMethod::quantile);

		normalization.set_no_series(params.mkmcParams.inputFilesPerSample.size());
		normalization.deserialize(StatisticsParams::NormalizationMethod::frequency_count, normalizationData);
		//norm.deserialize(StatisticsParams::NormalizationMethod::quantile, normalizationData);

		normalization.initialize();

		std::string kmerSequence;
		std::vector<size_t> matrixEntry;
		std::vector<double> normEntry;
		matrixEntry.resize(params.mkmcParams.inputFilesPerSample.size());
		normEntry.resize(params.mkmcParams.inputFilesPerSample.size());

		while (true)
		{
			if (!getLine(matrixFile, kmerSequence, matrixEntry))
				break;
			normalization.norm_entry(StatisticsParams::NormalizationMethod::frequency_count, matrixEntry, normEntry);
			putLine(frequenceNormFile, kmerSequence, normEntry);

			double pearson = refresh::correlation::pearson(normEntry.begin(), normEntry.end(), phenotype.begin());
			putLine(pearsonFile, kmerSequence, { pearson });

			refresh::correlation corr;

			double spearman = corr.spearman(normEntry.begin(), normEntry.end(), phenotype.begin());
			putLine(spearmanFile, kmerSequence, { spearman });

			double kendall = corr.kendall_tau(normEntry.begin(), normEntry.end(), phenotype.begin());
			putLine(kendallFile, kmerSequence, { kendall });
		}
	}
}


void StatisticsGenerator::generateStatisticsParallel()
{
	tasksData.reserve(params.stage1Params.GetNBins());
	for (uint32_t i = 0; i < params.stage1Params.GetNBins(); ++i)
	{
		tasksData.push_back(TaskData{ i });
	}

	readNormalizationDump(normalizationData, params.statisticsParams.normFrequencyFileTmp);
	readPhenotype(phenotype);

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
