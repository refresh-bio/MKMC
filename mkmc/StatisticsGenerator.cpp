#include "StatisticsGenerator.h"
#include "MatrixStats.h"
#include <algorithm>



void StatisticsGenerator::fillTaskData()
{
	try
	{
		matrixMetadataReader = std::make_unique<kmcdb::MetadataReader>(params.mkmcParams.outputFilesTemplate + ".kmcdb", false);
		matrixReader = std::make_unique<kmcdb::ReaderSortedPlainForListing<uint64_t>>(*matrixMetadataReader);

		kmcdb::Config config;
		config.num_bins = matrixMetadataReader->GetConfig().num_bins;
		config.signature_len = matrixMetadataReader->GetConfig().signature_len;
		config.signature_selection_scheme = matrixMetadataReader->GetConfig().signature_selection_scheme;
		config.signature_to_bin_mapping = matrixMetadataReader->GetConfig().signature_to_bin_mapping;
		config.kmer_len = matrixMetadataReader->GetConfig().kmer_len;
		config.num_samples = matrixMetadataReader->GetConfig().num_samples;
		config.num_bytes_single_value = { sizeof(out_kmcdb_value_type) };

		config.num_samples += params.statisticsParams.correlationMethods.size(); //I will add this correlations as a new columns
		//this make sense because those all of the same type (currently double), mkokot_TODO: remember to set appropriate col names!

		kmcdb::ConfigSortedPlain representation_config{};

		std::vector<std::string> sample_names;
		matrixReader->GetSampleNames(sample_names);
		assert(!sample_names.empty());

		auto is_method = [&](StatisticsParams::CorrelationMethod method)
			{
				const auto& corMeths = params.statisticsParams.correlationMethods;
				return std::find(corMeths.begin(), corMeths.end(), method) != corMeths.end();
			};

		if (is_method(StatisticsParams::CorrelationMethod::Pearson))
			sample_names.emplace_back("pearson_cor");

		if (is_method(StatisticsParams::CorrelationMethod::Spearman))
			sample_names.emplace_back("spearman_cor");

		if (is_method(StatisticsParams::CorrelationMethod::Kendall))
			sample_names.emplace_back("kendall_cor");

		kmcdbWriter = std::make_unique<kmcdb::WriterSortedPlain<double>>(
			config,
			representation_config,
			params.mkmcParams.outputFilesTemplate + "_norm+cor.kmcdb",
			params.mkmcParams.outputFilesTemplate + ".kmcdb",
			sample_names);
	}
	catch (const std::runtime_error& ex)
	{
		std::cerr << "Error: " << ex.what() << std::endl;
		exit(1);
	}
	tasksData.reserve(params.stage1Params.GetNBins());
	std::vector<uint64_t> nOutputKmersPerBin;
	nOutputKmersPerBin.reserve(params.stage1Params.GetNBins());
	
	for (uint32_t i = 0; i < params.stage1Params.GetNBins(); ++i)
	{
		tasksData.push_back(TaskData{ i });
		nOutputKmersPerBin.push_back(matrixReader->GetBin(i)->GetBinMetadata().total_kmers);
	}

	progress_bar = std::make_unique<ProgressBar>(
		params.mkmcParams.verbosity_level == 0 ? 0 : std::accumulate(nOutputKmersPerBin.begin(), nOutputKmersPerBin.end(), 0ull),
		"Computing statistics",
		std::cerr,
		params.mkmcParams.verbosity_level == 0);

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
		auto out_bin = kmcdbWriter->GetBin(taskData.binId);

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
		std::ofstream pearsonFile, spearmanFile, kendallFile;

		refresh::normalization_work<uint64_t, double> normalization;
		normalization.register_method(params.statisticsParams.normalizationMethod);
		normalization.set_no_series(params.mkmcParams.samples.size());
		normalization.deserialize(params.statisticsParams.normalizationMethod, normalizationData);

		normalization.initialize();

		refresh::correlation correlation;

		std::vector<uint64_t> matrixEntry;
		std::vector<double> outEntry;
		std::ptrdiff_t num_samples = static_cast<std::ptrdiff_t>(params.mkmcParams.samples.size());
		matrixEntry.resize(num_samples);
		outEntry.resize(num_samples);

		ProgressBarUpdater progress_bar_updater(*progress_bar, (std::max)(1ull, progress_bar->GetTotal() / 100ull));

		auto kmer_len = params.stage1Params.GetKmerLen();
		std::string kmerSequence(kmer_len, ' ');

		kmcdb::DispatchKmerSize<MAX_K>(kmer_len, [&](auto SIZE) {
			kmcdb::CKmer<SIZE> kmer;
			while (bin->NextKmer(kmer, matrixEntry.data()))
			{
				kmer.to_string(kmer_len, kmerSequence.data());
				normalization.norm_entry(params.statisticsParams.normalizationMethod, matrixEntry, outEntry);

				if (generatePearson)
				{
					const double pearson = refresh::correlation::pearson(
						outEntry.begin(),
						outEntry.begin() + num_samples,
						phenotype.begin());

					outEntry.push_back(pearson);
				}
				if (generateSpearman)
				{
					const double spearman = correlation.spearman(
						outEntry.begin(),
						outEntry.begin() + num_samples,
						phenotype.begin());

					outEntry.push_back(spearman);
				}
				if (generateKendall)
				{
					const double kendall = refresh::correlation::kendall_tau(
						outEntry.begin(),
						outEntry.begin() + num_samples,
						phenotype.begin());

					outEntry.push_back(kendall);
				}

				++progress_bar_updater;

				out_bin->AddKmer(kmer, outEntry.data());
			}
		});
	}
}


void StatisticsGenerator::generateStatisticsParallel()
{
	fillTaskData();

	if (params.statisticsParams.generateNormalization)
	{
		MatrixStatsReader stats_reader(params.mkmcParams.outputFilesTemplate + ".stats");
		bool success = false;
		if (params.statisticsParams.normalizationMethod == StatisticsParams::NormalizationMethod::frequency_count)
			success = stats_reader.Get(params.statisticsParams.normFrequencyStreamName, normalizationData);
		else if (params.statisticsParams.normalizationMethod == StatisticsParams::NormalizationMethod::quantile)
			success = stats_reader.Get(params.statisticsParams.normQuantileStreamName, normalizationData);

		if (!success)
		{
			std::cerr << "Error: cannot read normalization data\n";
			exit(1);
		}
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
