#include "StatisticsGenerator.h"
#include "MatrixStats.h"
#include "DumpWriter.h"
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
		config.num_samples += params.statisticsParams.generateEntropy ? 1 : 0;
		config.num_samples += params.statisticsParams.classificationMethods.size();
		//this make sense because those all of the same type (currently double)

		kmcdb::ConfigSortedPlain representation_config{};

		std::vector<std::string> sample_names;
		matrixReader->GetSampleNames(sample_names);
		assert(!sample_names.empty());

		if (statisticsToGeneration.pearson)
			sample_names.emplace_back("pearson_cor");

		if (statisticsToGeneration.spearman)
			sample_names.emplace_back("spearman_cor");

		if (statisticsToGeneration.kendall)
			sample_names.emplace_back("kendall_cor");

		if (statisticsToGeneration.entropy)
			sample_names.emplace_back("entropy");

		if (statisticsToGeneration.tTest)
			sample_names.emplace_back("ttest_analysis");

		if (statisticsToGeneration.snr)
			sample_names.emplace_back("snr_analysis");

		if (statisticsToGeneration.wilcoxonRankSum)
			sample_names.emplace_back("wrs_analysis");

		if (statisticsToGeneration.dids)
			sample_names.emplace_back("dids_analysis");

		if (statisticsToGeneration.anova)
			sample_names.emplace_back("anova_analysis");

		kmcdbWriter = std::make_unique<kmcdb::WriterSortedPlain<double>>(
			config,
			representation_config,
			params.mkmcParams.outputStatsBinFile,
			params.mkmcParams.outputBinFile,
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



StatisticsGenerator::StatisticsGenerator(Params& params) :
	params(params),
	tasksPool(tasksData),
	correlationPhenotype(params.phenotypes.correlationPhenotype.getPhenotype()),
	differentialAnalysisPhenotype(params.phenotypes.differentialAnalysisPhenotype.getMappedPhenotype()),
	differentialAnalysisNClasses(params.phenotypes.differentialAnalysisPhenotype.getClassesNumber())
{

	auto is_correlation_method = [&](StatisticsParams::CorrelationMethod method)
	{
		const auto& corMeths = params.statisticsParams.correlationMethods;
		return std::find(corMeths.begin(), corMeths.end(), method) != corMeths.end();
	};

	auto is_differential_analysis_method = [&](StatisticsParams::DifferentialAnalysisMethod method)
	{
		const auto& analysisMeths = params.statisticsParams.classificationMethods;
		return std::find(analysisMeths.begin(), analysisMeths.end(), method) != analysisMeths.end();
	};

	statisticsToGeneration.pearson = is_correlation_method(StatisticsParams::CorrelationMethod::Pearson);
	statisticsToGeneration.spearman = is_correlation_method(StatisticsParams::CorrelationMethod::Spearman);
	statisticsToGeneration.kendall = is_correlation_method(StatisticsParams::CorrelationMethod::Kendall);

	statisticsToGeneration.entropy = params.statisticsParams.generateEntropy;

	statisticsToGeneration.tTest = is_differential_analysis_method(StatisticsParams::DifferentialAnalysisMethod::TTest);
	statisticsToGeneration.snr = is_differential_analysis_method(StatisticsParams::DifferentialAnalysisMethod::SNR);
	statisticsToGeneration.wilcoxonRankSum = is_differential_analysis_method(StatisticsParams::DifferentialAnalysisMethod::WilcoxonRankSum);

	statisticsToGeneration.dids = is_differential_analysis_method(StatisticsParams::DifferentialAnalysisMethod::DIDS);
	statisticsToGeneration.anova = is_differential_analysis_method(StatisticsParams::DifferentialAnalysisMethod::ANOVA);

	statisticsToGeneration.differentialAnalysis = statisticsToGeneration.tTest || statisticsToGeneration.snr || statisticsToGeneration.wilcoxonRankSum || statisticsToGeneration.dids || statisticsToGeneration.anova;
}



void StatisticsGenerator::operator()()
{
	TaskData taskData;
	while (tasksPool.getTask(taskData))
	{
		auto bin = matrixReader->GetBin(taskData.binId);
		auto out_bin = kmcdbWriter->GetBin(taskData.binId);

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

		std::vector<uint64_t> matrixEntry;
		std::vector<double> outEntry;
		std::ptrdiff_t num_samples = static_cast<std::ptrdiff_t>(params.mkmcParams.samples.size());
		matrixEntry.resize(num_samples);
		outEntry.resize(num_samples);

		const uint32_t kmerLength = params.stage1Params.GetKmerLen();

		std::unique_ptr<OutputBuffer> pearsonOutputBuffer, spearmanOutputBuffer, kendallOutputBuffer;
		std::unique_ptr<OutputBuffer> entropyOutputBuffer;
		std::unique_ptr<OutputBuffer> tTestOutputBuffer, snrOutputBuffer, wilcoxonRankSumOutputBuffer;
		std::unique_ptr<OutputBuffer> didsOutputBuffer, anovaOutputBuffer;
		if (statisticsToGeneration.pearson)
			pearsonOutputBuffer = std::make_unique<OutputBuffer>(*writers.pearson, kmerLength);
		if (statisticsToGeneration.spearman)
			spearmanOutputBuffer = std::make_unique<OutputBuffer>(*writers.spearman, kmerLength);
		if (statisticsToGeneration.kendall)
			kendallOutputBuffer = std::make_unique<OutputBuffer>(*writers.kendall, kmerLength);

		if (statisticsToGeneration.entropy)
			entropyOutputBuffer = std::make_unique<OutputBuffer>(*writers.entropy, kmerLength);

		if (statisticsToGeneration.tTest)
			tTestOutputBuffer = std::make_unique<OutputBuffer>(*writers.tTest, kmerLength);
		if (statisticsToGeneration.snr)
			snrOutputBuffer = std::make_unique<OutputBuffer>(*writers.snr, kmerLength);
		if (statisticsToGeneration.wilcoxonRankSum)
			wilcoxonRankSumOutputBuffer = std::make_unique<OutputBuffer>(*writers.wilcoxonRankSum, kmerLength);

		if (statisticsToGeneration.dids)
			didsOutputBuffer = std::make_unique<OutputBuffer>(*writers.dids, kmerLength);
		if (statisticsToGeneration.anova)
			anovaOutputBuffer = std::make_unique<OutputBuffer>(*writers.anova, kmerLength);
		
		auto storeMethod = []<unsigned SIZE, typename VALUE_T>(const kmcdb::CKmer<SIZE>&kmer, uint64_t kmer_len, const VALUE_T cnt, char* out) -> size_t
		{
			kmer.to_string(kmer_len, out, '\t');
			size_t res = kmer_len + 1;
			out += kmer_len + 1;

			auto store_single_value = [&](const VALUE_T& val, char term)
			{
				size_t r{};
				if constexpr (std::is_integral_v<VALUE_T>)
					r = refresh::int_to_pchar(val, out, term);
				else if constexpr (std::is_floating_point_v<VALUE_T>)
				{
					if (std::isnan(val))
					{
						out[0] = 'n';
						out[1] = 'a';
						out[2] = 'n';
						out[3] = term;
						r = 4;
					}
					else if (std::isinf(val))
					{
						if (val < static_cast<VALUE_T>(0)) {
							out[0] = '-';
							out[1] = 'i';
							out[2] = 'n';
							out[3] = 'f';
							out[4] = term;
							r = 5;
						}
						else
						{
							out[0] = 'i';
							out[1] = 'n';
							out[2] = 'f';
							out[3] = term;
							r = 4;
						}
					}
					else
						r = refresh::real_to_pchar(val, out, 6, term);

				}
				else
				{
					static_assert(!sizeof(VALUE_T), "Unsupported type");
				}
				out += r;
				res += r;
			};

			store_single_value(cnt, '\n');

			return res;
		};

		ProgressBarUpdater progress_bar_updater(*progress_bar, (std::max)(1ull, progress_bar->GetTotal() / 100ull));

		auto kmer_len = params.stage1Params.GetKmerLen();
		std::string kmerSequence(kmer_len, ' ');

		kmcdb::DispatchKmerSize<MAX_K>(kmer_len, [&](auto SIZE) {
			kmcdb::CKmer<SIZE> kmer;
			while (bin->NextKmer(kmer, matrixEntry.data()))
			{
				kmer.to_string(kmer_len, kmerSequence.data());
				normalization.norm_entry(params.statisticsParams.normalizationMethod, matrixEntry, outEntry);

				if (statisticsToGeneration.pearson)
				{
					const double pearson = refresh::correlation::pearson_n(
						outEntry.begin(),
						correlationPhenotype.begin(),
						num_samples);

					outEntry.push_back(pearson);
					pearsonOutputBuffer->StoreKmer(kmer, kmerLength, pearson, storeMethod);
				}
				if (statisticsToGeneration.spearman)
				{
					const double spearman = correlation.spearman_n(
						outEntry.begin(),
						correlationPhenotype.begin(),
						num_samples);

					outEntry.push_back(spearman);
					spearmanOutputBuffer->StoreKmer(kmer, kmerLength, spearman, storeMethod);
				}
				if (statisticsToGeneration.kendall)
				{
					const double kendall = refresh::correlation::kendall_tau_n(
						outEntry.begin(),
						correlationPhenotype.begin(),
						num_samples);

					outEntry.push_back(kendall);
					kendallOutputBuffer->StoreKmer(kmer, kmerLength, kendall, storeMethod);
				}

				if (statisticsToGeneration.entropy)
				{
					const double entropy = entropyObj.entropy_n(
						outEntry.begin(),
						num_samples);

					outEntry.push_back(entropy);
					entropyOutputBuffer->StoreKmer(kmer, kmerLength, entropy, storeMethod);
				}
				if (statisticsToGeneration.differentialAnalysis)
				{
					if (statisticsToGeneration.tTest)
					{
						const double tTestPValue = statistics.t_test_n(
							outEntry.begin(),
							differentialAnalysisPhenotype.begin(),
							num_samples).p_value;

						outEntry.push_back(tTestPValue);
						tTestOutputBuffer->StoreKmer(kmer, kmerLength, tTestPValue, storeMethod);
					}
					if (statisticsToGeneration.snr)
					{
						const double snr = statistics.SNR_test_n(
							outEntry.begin(),
							differentialAnalysisPhenotype.begin(),
							num_samples);

						outEntry.push_back(snr);
						snrOutputBuffer->StoreKmer(kmer, kmerLength, snr, storeMethod);
					}
					if (statisticsToGeneration.wilcoxonRankSum)
					{
						const double wilcoxonRankSumPValue = statistics.mann_whitney_U_test_n(
							outEntry.begin(),
							differentialAnalysisPhenotype.begin(),
							num_samples).p_value;

						outEntry.push_back(wilcoxonRankSumPValue);
						wilcoxonRankSumOutputBuffer->StoreKmer(kmer, kmerLength, wilcoxonRankSumPValue, storeMethod);
					}
					if (statisticsToGeneration.dids)
					{
						const double dids = scorer.dids_n(
							outEntry.begin(),
							differentialAnalysisPhenotype.begin(),
							differentialAnalysisNClasses,
							num_samples);

						outEntry.push_back(dids);
						didsOutputBuffer->StoreKmer(kmer, kmerLength, dids, storeMethod);
					}
					if (statisticsToGeneration.anova)
					{
						const double anovaPValue = scorer.anova_n(
							outEntry.begin(),
							differentialAnalysisPhenotype.begin(),
							differentialAnalysisNClasses,
							num_samples).p_value;

						outEntry.push_back(anovaPValue);
						anovaOutputBuffer->StoreKmer(kmer, kmerLength, anovaPValue, storeMethod);
					}
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
		MatrixStatsReader stats_reader(params.mkmcParams.normStatsBinFile);
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


	const bool multiThreadedGeneration = params.mkmcParams.nThreads > 1;

	if (statisticsToGeneration.pearson)
	{
		writers.pearson = std::make_unique<DumpWriter>(params.mkmcParams.outputFilePearson, multiThreadedGeneration);
		writers.pearson->StoreHeader({ "pearson_cor" });
	}
	if (statisticsToGeneration.spearman)
	{
		writers.spearman = std::make_unique<DumpWriter>(params.mkmcParams.outputFileSpearman, multiThreadedGeneration);
		writers.spearman->StoreHeader({ "spearman_cor" });
	}
	if (statisticsToGeneration.kendall)
	{
		writers.kendall = std::make_unique<DumpWriter>(params.mkmcParams.outputFileKendall, multiThreadedGeneration);
		writers.kendall->StoreHeader({ "kendall_cor" });
	}
	if (statisticsToGeneration.entropy)
	{
		writers.entropy = std::make_unique<DumpWriter>(params.mkmcParams.outputFileEntropy, multiThreadedGeneration);
		writers.entropy->StoreHeader({ "entropy" });
	}
	if (statisticsToGeneration.tTest)
	{
		writers.tTest = std::make_unique<DumpWriter>(params.mkmcParams.outputFileTTest, multiThreadedGeneration);
		writers.tTest->StoreHeader({ "ttest_analysis_p_val" });
	}
	if (statisticsToGeneration.snr)
	{
		writers.snr = std::make_unique<DumpWriter>(params.mkmcParams.outputFileSNR, multiThreadedGeneration);
		writers.snr->StoreHeader({ "snr_analysis" });
	}
	if (statisticsToGeneration.wilcoxonRankSum)
	{
		writers.wilcoxonRankSum = std::make_unique<DumpWriter>(params.mkmcParams.outputFileWilcoxonRankSum, multiThreadedGeneration);
		writers.wilcoxonRankSum->StoreHeader({ "wrs_analysis_p_val" });
	}
	if (statisticsToGeneration.dids)
	{
		writers.dids = std::make_unique<DumpWriter>(params.mkmcParams.outputFileDIDS, multiThreadedGeneration);
		writers.dids->StoreHeader({ "dids_analysis" });
	}
	if (statisticsToGeneration.anova)
	{
		writers.anova = std::make_unique<DumpWriter>(params.mkmcParams.outputFileANOVA, multiThreadedGeneration);
		writers.anova->StoreHeader({ "anova_analysis_p_val" });
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
