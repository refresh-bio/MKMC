#include "StatisticsGenerator.h"
#include <algorithm>



void StatisticsGenerator::fillTaskData()
{
	try
	{
		matrixMetadataReader = std::make_unique<kmcdb::MetadataReader>(params.mkmcParams.outputFilesTemplate, false);
		matrixReader = std::make_unique<kmcdb::ReaderSortedPlainForListing<uint64_t>>(*matrixMetadataReader);
	}
	catch (const std::runtime_error& ex)
	{
		std::cerr << "Error: " << ex.what() << std::endl;
		exit(1);
	}
	tasksData.reserve(params.stage1Params.GetNBins());
	for (uint32_t i = 0; i < params.stage1Params.GetNBins(); ++i)
	{
		tasksData.push_back(TaskData{ i });
	}
	std::vector<uint64_t> nOutputKmersPerBin;
	readDump(nOutputKmersPerBin, params.statisticsParams.statsNOutputKmers);
	std::sort(tasksData.begin(), tasksData.end(), [&](const TaskData& a, const TaskData& b) { return nOutputKmersPerBin[a.binId] > nOutputKmersPerBin[b.binId]; });
}



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

	if (phenotype.size() != params.mkmcParams.samples.size())
	{
		std::cerr << "Error: a phenotype size in a file " << params.statisticsParams.phenotypeFile  <<  " (" << phenotype.size() << ") is different than number of samples (" << params.mkmcParams.samples.size() << ")." << std::endl;
		exit(1);
	}
}



void StatisticsGenerator::operator()()
{
	TaskData taskData;
	while (tasksPool.getTask(taskData))
	{
		auto bin = matrixReader->GetBin(taskData.binId);

		std::ofstream normFile(params.mkmcParams.outputFilesNorm[taskData.binId]);
		if (!normFile.is_open())
		{
			std::cerr << "Error: cannot open " << params.mkmcParams.outputFilesNorm[taskData.binId] << "." << std::endl;
			exit(1);
		}

		bool generatePearson = false, generateSpearman = false, generateKendall = false;
		for (auto method : params.statisticsParams.correlationMethods)
		{
			if (method == StatisticsParams::CorrelationMethod::Pearson)
				generatePearson = true;
			else if (method == StatisticsParams::CorrelationMethod::Spearman)
				generateSpearman = true;
			else if (method == StatisticsParams::CorrelationMethod::Kendall)
				generateKendall = true;
		}

		//mkokot_TODO: ok, for now I will just generate text file, but later i will write to a common kmcdb
		//instead of norm file etc, so sample names will be taken directly from input kmcdb
		normFile << "k-mer\t";
		for (const Sample& sample : params.mkmcParams.samples)
			normFile << sample.name << '\t';

		normFile << '\n';

		std::ofstream pearsonFile, spearmanFile, kendallFile;
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

		refresh::normalization_work<uint64_t, double> normalization;
		normalization.register_method(params.statisticsParams.normalizationMethod);
		normalization.set_no_series(params.mkmcParams.samples.size());
		normalization.deserialize(params.statisticsParams.normalizationMethod, normalizationData);

		normalization.initialize();

		refresh::correlation correlation;

		std::vector<uint64_t> matrixEntry;
		std::vector<double> normEntry;
		matrixEntry.resize(params.mkmcParams.samples.size());
		normEntry.resize(params.mkmcParams.samples.size());

		ProgressBarUpdater progress_bar_updater(progress_bar, (std::max)(1ull, totAllKmers / 100ull));

		auto kmer_len = params.stage1Params.GetKmerLen();
		std::string kmerSequence(kmer_len, ' ');

		kmcdb::DispatchKmerSize<MAX_K>(kmer_len, [&](auto SIZE) {
			kmcdb::CKmer<SIZE> kmer;
			while (bin->NextKmer(kmer, matrixEntry.data()))
			{
				kmer.to_string(kmer_len, kmerSequence.data());
				normalization.norm_entry(params.statisticsParams.normalizationMethod, matrixEntry, normEntry);
				putLine(normFile, kmerSequence, normEntry);

				if (generatePearson)
				{
					const double pearson = refresh::correlation::pearson(normEntry.begin(), normEntry.end(), phenotype.begin());
					putLine(pearsonFile, kmerSequence, { pearson });
				}
				if (generateSpearman)
				{
					const double spearman = correlation.spearman(normEntry.begin(), normEntry.end(), phenotype.begin());
					putLine(spearmanFile, kmerSequence, { spearman });
				}
				if (generateKendall)
				{
					const double kendall = refresh::correlation::kendall_tau(normEntry.begin(), normEntry.end(), phenotype.begin());
					putLine(kendallFile, kmerSequence, { kendall });
				}

				++progress_bar_updater;
			}
		});
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

	if (!params.statisticsParams.correlationMethods.empty())
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
