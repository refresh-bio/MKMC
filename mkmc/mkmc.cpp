#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <map>
#include <functional>
#include "kmc_core/kmc_runner.h"
#include "parameters.h"
#include "KMCRunner.h"
#include "Merger.h"
#include "Finish.h"
#include "Start.h"
#include "Version.h"
#include "time.hpp"
#include "SequenceFilterInit.h"
#include "Logger.h"
#include "StatisticsGenerator.h"
#include "CLI11/include/CLI/CLI.hpp"



class MergerRunner
{
	Params& params;
	Timer& merger_timer;
public:
	MergerRunner(Params& params, Timer& merger_timer):
		params(params),
		merger_timer(merger_timer)
	{

	}
	template<unsigned SIZE>
	void Run()
	{
		std::vector<std::string> tasks = { params.mkmcParams.outputMatrixBinFile };
		const auto& outputFileTypes = params.mkmcParams.outputFileTypes;
		if (std::find(outputFileTypes.begin(), outputFileTypes.end(), OutputFileType::Matrix) != outputFileTypes.end())
			tasks.push_back(params.mkmcParams.outputMatrixFile);
		if (std::find(outputFileTypes.begin(), outputFileTypes.end(), OutputFileType::FASTA) != outputFileTypes.end())
			tasks.push_back(params.mkmcParams.outputFASTAFile);

		std::cerr << "Starting merging samples and dumping to " << (tasks.size() > 1 ? "files " : "file ") << MessagesUtilities::generateStartingSentence(tasks) << "..." << std::endl << std::endl;

		Merger<SIZE> merger(params);
		merger_timer.startTimer();
		merger.mergeParallel();
		merger_timer.stopTimer();
	}
};



void configureArguments(int argc, char** argv, Params& params, CLI::App& app)
{
	DefaultKMCParams& defaultKMCParams = params.defaultKMCParams;
	KMC::Stage1Params& stage1Params = params.stage1Params;
	MKMCParams& mkmcParams = params.mkmcParams;
	FilterParams& filterParams = params.filterParams;
	StatisticsParams& statisticsParams = params.statisticsParams;
	Phenotypes& phenotypes = params.phenotypes;

	// set KMC defaults
	stage1Params.SetKmerLen(defaultKMCParams.k);
	stage1Params.SetCutoffMin(defaultKMCParams.ci);
	stage1Params.SetCutoffMax(defaultKMCParams.cx);
	stage1Params.SetCounterMax(defaultKMCParams.cs);

	CLI::Option* p = nullptr, * n = nullptr, * cor = nullptr, * differentialAnalysis = nullptr, *pvalCorr, * c = nullptr;

	app.add_option("input_samples_file", mkmcParams.inputFileName, "file with a list of samples names with input files names in specified (-f parameter) format (gzipped or not)")->required()->check(CLI::ExistingFile);
	app.add_option("output_files_template", mkmcParams.outputFilesTemplate, "template (prefix) of output files names")->required();
	app.add_option("temp_dir", mkmcParams.tmpPath, "a directory where temporary files will be stored")->required();

	std::function<void(const uint32_t&)> kCallback = [&](const uint32_t& k)
	{
		stage1Params.SetKmerLen(k);
	};
	app.add_option_function("-k", kCallback, "k-mer length")->check(CLI::Range(KMC::CfgConsts::min_k, KMC::CfgConsts::max_k))->default_val(defaultKMCParams.k);

	app.add_flag("--tot_cnt", mkmcParams.totCntGeneration, "generate samples counts sums file");

	CLI::Option_group* filteringGroup = app.add_option_group("k-mers filtering");
	filteringGroup->add_option("--thr", filterParams.minCountThreshold, "filter out k-mers occuring less than specified number of times...")->check(CLI::PositiveNumber)->default_val(filterParams.minCountThreshold);
	filteringGroup->add_option("--thr_rat", filterParams.minKmersAboveThresholdRatio, "... in a specified ratio of the input files (see example)")->check(CLI::Range(0.0, 1.0))->default_val(filterParams.minKmersAboveThresholdRatio);

	std::function<void(const decltype(filterParams.inputKmersSequencesToFilterOut)&)> fltCallback = [&](const decltype(filterParams.inputKmersSequencesToFilterOut)& fileName)
	{
		filterParams.inputKmersSequencesToFilterOut = fileName;
		filterParams.filterKmersSequences = true;
	};
	filteringGroup->add_option_function("--flt", fltCallback, "keep k-mers present in a specified file (FASTA or a set of the k-mers, one in each line) only; if -b is not set, the k-mers are converted to canonical form")->check(CLI::ExistingFile);

	CLI::Option_group* correlationGroup = app.add_option_group("correlation and normalization");

	std::map<std::string, StatisticsParams::NormalizationMethod> valuesMap{ {"deseq", StatisticsParams::NormalizationMethod::deseq2}, {"freq", StatisticsParams::NormalizationMethod::frequency_count }, {"q", StatisticsParams::NormalizationMethod::quantile } };
	std::function<void(const decltype(statisticsParams.normalizationMethod)&)> nCallback = [&](const decltype(statisticsParams.normalizationMethod)& normalizationMethod)
	{
		statisticsParams.normalizationMethod = normalizationMethod;
		statisticsParams.generateNormalization = true;
	};
	n = correlationGroup->add_option_function("-n", nCallback, "generate normalized counts (DESeq2/frequency count/quantile normalization)")->transform(CLI::CheckedTransformer(valuesMap, CLI::ignore_case));

	std::map<std::string, StatisticsParams::CorrelationMethod> correlationValuesMap{ {"pearson", StatisticsParams::CorrelationMethod::Pearson }, { "spearman", StatisticsParams::CorrelationMethod::Spearman }, {"kendall", StatisticsParams::CorrelationMethod::Kendall } };
	cor = correlationGroup->add_option("--cor", statisticsParams.correlationMethods, "compute correlation cofficients, basing on a phenotype file (Kendall Tau/Pearson/Spearman correlation)")->transform(CLI::CheckedTransformer(correlationValuesMap));

	std::function<void(const std::string&)> pCallback = [&](const std::string& fileName)
	{
		phenotypes.correlationPhenotype.setFileName(fileName);
	};
	p = correlationGroup->add_option_function("-p", pCallback, "set a phenotype file (a sequence of integers, one in each line)")->check(CLI::ExistingFile)->needs(cor);

	CLI::Option_group* diffGroup = app.add_option_group("differential k-mers analysis");
	typedef StatisticsParams::DifferentialAnalysisMethod DAMethod;
	std::map<std::string, DAMethod> differentialAnalysisValuesMap{ {"ttest", DAMethod::TTest }, {"snr", DAMethod::SNR }, {"wrs", DAMethod::WilcoxonRankSum }, {"dids", DAMethod::DIDS }, {"anova", DAMethod::ANOVA } };
	differentialAnalysis = diffGroup->add_option("--diff", statisticsParams.classificationMethods, "perform differential k-mers analysis (ANOVA, DIDS, Signal to Noise ratio, T-Test, Wilcoxon-rank sum (Mann-Whitney U test)); all except T-Test need -n; counts for T-Test are increased by 1 and logarithmized")->transform(CLI::CheckedTransformer(differentialAnalysisValuesMap));

	typedef StatisticsParams::DifferentialAnalysisCorrectionMethod CorrectionMethod;
	std::map<std::string, StatisticsParams::DifferentialAnalysisCorrectionMethod> differentialAnalysisCorrectionValuesMap{ { "b", CorrectionMethod::Bonferroni }, { "hb", CorrectionMethod::HolmBonferroni }, { "bh", CorrectionMethod::BenjaminiHochberg }, { "by", CorrectionMethod::BenjaminiYekutieli } };
	std::function<void(const decltype(statisticsParams.classificationPValueCorrection)&)> pcorrCallback = [&](const decltype(statisticsParams.classificationPValueCorrection)& classificationPValueCorrection)
	{
		statisticsParams.classificationPValueCorrection = classificationPValueCorrection;
		statisticsParams.correctPvalues = true;
	};
	pvalCorr = diffGroup->add_option_function("--pval_corr", pcorrCallback, "correct p-values of differential k-mers analysis (Bonferroni, Benjamini-Hochberg, Benjamini-Yekutieli, Holm-Bonferroni); store statistically significant k-mers also in separated files")->transform(CLI::CheckedTransformer(differentialAnalysisCorrectionValuesMap))->needs(differentialAnalysis);

	diffGroup->add_option("--max_corrected_pval", statisticsParams.maxCorrectedPval, "statistical significance for --pval_corr parameter")->check(CLI::Range(0.0, 1.0))->default_val(statisticsParams.maxCorrectedPval)->needs(pvalCorr);

	std::function<void(const std::string&)> cCallback = [&](const std::string& fileName)
	{
		phenotypes.differentialAnalysisPhenotype.setFileName(fileName);
	};
	c = diffGroup->add_option_function("-c", cCallback, "set a phenotype file for differential k-mers analysis (a sequence of natural numbers or text labels, one in each line)")->check(CLI::ExistingFile)->needs(differentialAnalysis);

	typedef StatisticsParams::DIDSMode DIDSMode;
	std::map<std::string, DIDSMode> didsModeValuesMap{ { "sqrt", DIDSMode::sqrt }, { "quadratic", DIDSMode::quadratic }, { "tanh", DIDSMode::tanh } };
	std::function<void(const decltype(statisticsParams.didsMode)&)> didsModeCallback = [&](const decltype(statisticsParams.didsMode)& didsMode)
	{
		statisticsParams.didsMode = didsMode;
		statisticsParams.didsModeUserDefined = true;
	};
	diffGroup->add_option_function("--dids-mode", didsModeCallback, "DIDS mode (x*x, square root, 1 + tanh(3x - 3)")->transform(CLI::CheckedTransformer(didsModeValuesMap))->default_val(statisticsParams.didsMode)->default_str("sqrt");


	CLI::Option_group* cvGroup = app.add_option_group("cross-validation");
	auto cv = cvGroup->add_flag("--cv", statisticsParams.cvParams.cv, "perform cross-validation for correlation")->needs(cor);
	cvGroup->add_option("--leave", statisticsParams.cvParams.p, "number of samples to leave in every test")->needs(cv)->default_val(statisticsParams.cvParams.p);
	std::function<void(const decltype(statisticsParams.cvParams.seed)&)> nCVSeed = [&](const decltype(statisticsParams.cvParams.seed)& seed)
	{
		statisticsParams.cvParams.seed = seed;
		statisticsParams.cvParams.seedUserDefined = true;
	};
	cvGroup->add_option_function("--cv-seed", nCVSeed, "random seed")->needs(cv)->default_val(statisticsParams.cvParams.seed);


	CLI::Option_group* otherStatsGroup = app.add_option_group("other statistical parameters");

	otherStatsGroup->add_flag("--entropy", statisticsParams.generateEntropy, "generate k-mers counts entropy; counts are increased by 1");

	std::function<void(const decltype(statisticsParams.nTop)&)> nTopCallback = [&](const decltype(statisticsParams.nTop)& nTop)
	{
		statisticsParams.nTop = nTop;
		statisticsParams.nTopUserDefined = true;
	};
	otherStatsGroup->add_option_function("--n_top", nTopCallback, "select a maximal number of top k-mers by statistics with no p-values (for correlations in terms of an absolute value) and store them in separate files; needs --cor or --diff")->default_val(statisticsParams.nTop);
	
	CLI::Option_group* dimReductionGroup = app.add_option_group("dimentionality reduction");

	auto umap = dimReductionGroup->add_flag("--umap", statisticsParams.runUMAP, "run dimentionality reduction on normalized matrix with UMAP")->needs(n);
	auto pca = dimReductionGroup->add_flag("--pca", statisticsParams.runPCA, "run dimentionality reduction on normalized matrix with PCA")->needs(n);

	std::function<void(const decltype(statisticsParams.nDimensionReduction)&)> dimensionsCallback = [&](const decltype(statisticsParams.nDimensionReduction)& dimensions)
	{
		statisticsParams.nDimensionReduction = dimensions;
		statisticsParams.nDimensionReductionUserDefined = true;
	};
	dimReductionGroup->add_option_function("--dimensions", dimensionsCallback, "number of output dimensions; needs --umap or --pca")->default_val(statisticsParams.nDimensionReduction);

	dimReductionGroup->add_option("--umap-local_connectivity", statisticsParams.umap_params.local_connectivity, "local_connectivity parameter")->needs(umap)->default_val(statisticsParams.umap_params.local_connectivity);
	dimReductionGroup->add_option("--umap-bandwidth", statisticsParams.umap_params.bandwidth, "bandwidth parameter")->needs(umap)->default_val(statisticsParams.umap_params.bandwidth);

	dimReductionGroup->add_option("--umap-mix_ratio", statisticsParams.umap_params.mix_ratio, "mix_ratio parameter")->needs(umap)->default_val(statisticsParams.umap_params.mix_ratio);
	dimReductionGroup->add_option("--umap-spread", statisticsParams.umap_params.spread, "spread parameter")->needs(umap)->default_val(statisticsParams.umap_params.spread);
	dimReductionGroup->add_option("--umap-min_dist", statisticsParams.umap_params.min_dist, "min_dist parameter")->needs(umap)->default_val(statisticsParams.umap_params.min_dist);
	dimReductionGroup->add_option("--umap-a", statisticsParams.umap_params.a, "a parameter")->needs(umap)->default_val(statisticsParams.umap_params.a);
	dimReductionGroup->add_option("--umap-b", statisticsParams.umap_params.b, "b parameter")->needs(umap)->default_val(statisticsParams.umap_params.b);
	dimReductionGroup->add_option("--umap-repulsion_strength", statisticsParams.umap_params.repulsion_strength, "repulsion_strength parameter")->needs(umap)->default_val(statisticsParams.umap_params.repulsion_strength);

	std::map<std::string, umappp::InitMethod> umapInitMethodValuesMap{ { "spectral", umappp::InitMethod::SPECTRAL }, { "spectral_only", umappp::InitMethod::SPECTRAL_ONLY }, { "random", umappp::InitMethod::RANDOM }, { "none", umappp::InitMethod::NONE } };
	dimReductionGroup->add_option("--umap-initialize", statisticsParams.umap_params.initialize, "initialize parameter")->transform(CLI::CheckedTransformer(umapInitMethodValuesMap))->needs(umap)->default_val(statisticsParams.umap_params.initialize)->default_str("spectral");

	dimReductionGroup->add_option("--umap-num_neighbors", statisticsParams.umap_params.num_neighbors, "num_neighbors parameter")->needs(umap)->default_val(statisticsParams.umap_params.num_neighbors);
	dimReductionGroup->add_option("--umap-num_epochs", statisticsParams.umap_params.num_epochs, "num_epochs parameter")->needs(umap)->default_val(statisticsParams.umap_params.num_epochs); // default -1
	dimReductionGroup->add_option("--umap-learning_rate", statisticsParams.umap_params.learning_rate, "learning_rate parameter")->needs(umap)->default_val(statisticsParams.umap_params.learning_rate);

	dimReductionGroup->add_option("--umap-negative_sample_rate", statisticsParams.umap_params.negative_sample_rate, "negative_sample_rate parameter")->needs(umap)->default_val(statisticsParams.umap_params.negative_sample_rate);
	dimReductionGroup->add_option("--umap-seed", statisticsParams.umap_params.seed, "seed parameter")->needs(umap)->default_val(statisticsParams.umap_params.seed);

	//this will be set with the "main" or "global" number of threads
	//umapGroup->add_option("--umap-num_threads", statisticsParams.umap_params.num_threads, "num_threads parameter of umap")->needs(umap)->default_val(statisticsParams.umap_params.num_threads);

	dimReductionGroup->add_option("--umap-parallel_optimization", statisticsParams.umap_params.parallel_optimization, "parallel_optimization parameter")->needs(umap)->default_val(statisticsParams.umap_params.parallel_optimization);

	std::map<std::string, refresh::pca<double>::computation_mode_t> pcaModeValuesMap{ { "svd", refresh::pca<double>::computation_mode_t::svd }, { "covariance", refresh::pca<double>::computation_mode_t::covariance } };
	dimReductionGroup->add_option("--pca-mode", statisticsParams.pca_mod, "PCA mode")->transform(CLI::CheckedTransformer(pcaModeValuesMap))->needs(pca)->default_val(statisticsParams.pca_mod)->default_str("svd");


	CLI::Option_group* optionalGroup = app.add_option_group("additional parameters");

	std::map<std::string, KMC::InputFileType> inputValuesMap{ {"fa", KMC::InputFileType::FASTA }, {"fq", KMC::InputFileType::FASTQ }, { "mf", KMC::InputFileType::MULTILINE_FASTA } };
	optionalGroup->add_option("-f", mkmcParams.inputFileType, "input format (FASTA, FASTQ or multi-FASTA); mixing files formats is not supported")->transform(CLI::CheckedTransformer(inputValuesMap, CLI::ignore_case))->default_val(mkmcParams.inputFileType)->default_str("fq");

	std::map<std::string, OutputFileType> outputValuesMap{ {"fa", OutputFileType::FASTA }, {"matrix", OutputFileType::Matrix } };
	optionalGroup->add_option("-o", mkmcParams.outputFileTypes, "output format (FASTA or matrix)")->transform(CLI::CheckedTransformer(outputValuesMap, CLI::ignore_case));

	std::function<void()> bCallback = [&]()
	{
		stage1Params.SetCanonicalKmers(false);
	};
	optionalGroup->add_flag_callback("-b", bCallback, "turn off transformation of k-mers into canonical form; applies both for input sequences and k-mers passed by --flt");

	std::function<void(const uint32_t&)> ciCallback = [&](const uint32_t& ci) // currently 32 bits
	{
		stage1Params.SetCutoffMin(static_cast<uint64_t>(ci));
	};
	optionalGroup->add_option_function("--ci", ciCallback, "exclude k-mers occurring less than specified number of times (if k-mer occurs less than --ci times in a sample, it gets counter 0, but for this sample only)")->check(CLI::PositiveNumber)->default_val(defaultKMCParams.ci);
	std::function<void(const uint32_t&)> cxCallback = [&](const uint32_t& cx) // currently 32 bits
	{
		stage1Params.SetCutoffMax(static_cast<uint64_t>(cx));
	};
	optionalGroup->add_option_function("--cx", cxCallback, "exclude counting k-mers occurring more than specified number of times (if k-mer occurs more than --cx times in a sample, it gets counter 0, but for this sample only)")->check(CLI::PositiveNumber)->default_val(static_cast<uint32_t>(defaultKMCParams.cx));

	std::function<void(const uint32_t&)> csCallback = [&](const uint32_t& cs) // currently 32 bits
	{
		stage1Params.SetCounterMax(static_cast<uint64_t>(cs));
	};
	optionalGroup->add_option_function("--cs", csCallback, "maximal value of a counter")->check(CLI::Range(2U, std::numeric_limits<uint32_t>::max()))->default_val(defaultKMCParams.cs);

	std::function<void(const decltype(mkmcParams.nKMCWorkers)&)> wrkCallback = [&](const decltype(mkmcParams.nKMCWorkers)& nKMCWorkers)
	{
		mkmcParams.nKMCWorkers = nKMCWorkers;
		mkmcParams.nKMCWorkersUserSet = true;
	};
	optionalGroup->add_option_function("--wrk", wrkCallback, "number of parallel k-mer counting tasks")->default_val(mkmcParams.nKMCWorkers);

	optionalGroup->add_option("-t", mkmcParams.nThreads, "number of threads")->default_val(mkmcParams.nThreads);

	std::function<void(const uint32_t& maxRamGB)> mCallback = [&](const uint32_t& maxRamGB)
	{
		mkmcParams.maxRamGB = maxRamGB;
		mkmcParams.maxRamGBUserDefined = true;
	};
	optionalGroup->add_option_function("-m", mCallback, "max amount of RAM in GB; practically works only if -r is not set")->check(CLI::Range(2, 1024))->default_val(mkmcParams.maxRamGB);

	std::function<void()> rCallback = [&]()
	{
		stage1Params.SetRamOnlyMode(true);
	};
	optionalGroup->add_flag_callback("-r", rCallback, "RAM only mode for k-mer counting");

	std::function<void()> vCallback = [&]()
	{
		mkmcParams.verbosity_level++;
	};
	optionalGroup->add_flag_callback("-v", vCallback, "verbose mode, shows progress");

	CLI::Option_group* debugGroup = app.add_option_group("debug parameters");
	debugGroup->add_flag("--keep", mkmcParams.keepTmpFiles, "keep temporary files and binary results file");

	debugGroup->add_option("--on", mkmcParams.nKMCBins, "number of internal bins, modify carefully")->check(CLI::PositiveNumber)->default_val(mkmcParams.nKMCBins);

	debugGroup->add_flag("--generate_snr_for_unnormalized_data", mkmcParams.generateForNonNormalized, "generate Signal to Noise ratio also for unnormalized counts");

	cor->needs(n)->needs(p);
	differentialAnalysis->needs(c);

	app.footer("Warning: k-mers order in output files is not specified and may vary between runnings.\n\n"
		"Example: to run MKMC, type:\n"
		"    ./mkmc -k 20 --thr_rat 0.5 input_files_list.txt output tmp\n"
		"It will generate a matrix of 20-mers occurring in at least a half of the input files.\n"
		"    ./mkmc -k 20 --thr 2 --thr_rat 0.5 input_files_list.txt output tmp\n"
		"It will generate a matrix of 20-mers occurring at least twice in at least a half of the input files.\n\n"
		"input_files_list.txt example:\n"
		"    killifishretina1 kfA_1.fastq.gz kfA_2.fastq.gz\n"
		"    killifishretina2 kfB.fastq.gz");
}



bool checkAndPrintArgumentsErrors(const Params& params)
{
	const StatisticsParams statisticsParams = params.statisticsParams;

	auto isDAMethod = [&](StatisticsParams::DifferentialAnalysisMethod method)
	{
		return std::find(statisticsParams.classificationMethods.begin(), statisticsParams.classificationMethods.end(), method) != statisticsParams.classificationMethods.end();
	};


	if (!statisticsParams.generateNormalization &&
		(statisticsParams.classificationMethods.size() > 1 || statisticsParams.classificationMethods.size() == 1 && statisticsParams.classificationMethods.front() != StatisticsParams::DifferentialAnalysisMethod::TTest))
	{
		std::cerr << "Error: Differential analysis methods (except T-Test) require normalization (-n)\n";
		return true;
	}

	if (statisticsParams.nDimensionReductionUserDefined &&
		!(statisticsParams.runPCA || statisticsParams.runUMAP))
	{
		std::cerr << "Error: Number of dimensions (--dimensions) requires dimensionality reduction algorithm (--umap or --pca)\n";
		return true;
	}

	if (statisticsParams.correctPvalues &&
		!isDAMethod(StatisticsParams::DifferentialAnalysisMethod::ANOVA) &&
		!isDAMethod(StatisticsParams::DifferentialAnalysisMethod::TTest) &&
		!isDAMethod(StatisticsParams::DifferentialAnalysisMethod::WilcoxonRankSum))
	{
		std::cerr << "Error: --pval_corr requires differential k-mers analysis with ANOVA, T-Test, or Wilcoxon-rank sum (Mann-Whitney U test) (--diff)";
		return true;
	}

	if (statisticsParams.didsModeUserDefined && !isDAMethod(StatisticsParams::DifferentialAnalysisMethod::DIDS))
	{
		std::cerr << "Error: --dids-mode requires differential k-mers analysis with DIDS (--diff)";
		return true;
	}

	if (statisticsParams.nTopUserDefined &&
		statisticsParams.correlationMethods.empty() &&
		!statisticsParams.generateEntropy &&
		!isDAMethod(StatisticsParams::DifferentialAnalysisMethod::SNR) &&
		!isDAMethod(StatisticsParams::DifferentialAnalysisMethod::DIDS))
	{
		std::cerr << "Error: --n_top requires correlation (--cor) or entropy (--entropy) or differential k-mers analysis with SNR or DIDS (--diff)";
		return true;
	}

	if (statisticsParams.nTopUserDefined &&
		statisticsParams.nTop < 1)
	{
		std::cerr << "Error: --n_top has to be at least 1";
		return true;
	}

	return false;
}



bool checkAndPrintParamsFromFilesErrors(const Params& params)
{
	if ((params.statisticsParams.runPCA || params.statisticsParams.runUMAP) &&
		(params.statisticsParams.nDimensionReduction < 1 || params.statisticsParams.nDimensionReduction >= params.mkmcParams.samples.size()))
	{
		std::cerr << "Error: Number of dimensions (--dimensions) must be at least 1 and lower than number of samples\n";
		return true;
	}

	if (params.statisticsParams.cvParams.p == 0 || params.statisticsParams.cvParams.p >= params.mkmcParams.samples.size() || params.mkmcParams.samples.size() % params.statisticsParams.cvParams.p != 0)
	{
		std::cerr << "Error: Number of samples to leave in cross-validation (--leave) has to be positive and be a factor of a number of samples.\n";
		return true;
	}
	return false;
}


bool checkAndPrintArgumentsWarnings(const Params& params)
{
	bool result = false;
	if (params.mkmcParams.samples.size() <= 8 && std::find(params.statisticsParams.classificationMethods.begin(), params.statisticsParams.classificationMethods.end(), StatisticsParams::DifferentialAnalysisMethod::WilcoxonRankSum) != params.statisticsParams.classificationMethods.end())
	{
		std::cerr << "Warning: Wilcoxon-rank sum (Mann-Whitney U test) uses approximate algorithm, thus for less than 9 samples its results may be slightly different than in e.g. SciPy.\n";
		result = true;
	}

	if (params.statisticsParams.cvParams.seedUserDefined && params.statisticsParams.cvParams.p == 1)
	{
		std::cerr << "Warning: as --leave parameter is set to 1, LOOCV will be performed, which does not need randomness (--cv-seed parameter will be ignored).\n";
		result = true;
	}

	return result;
}



//----------------------------------------------------------------------------------
// Main function
int main(int argc, char** argv)
{
	Params params;

	// CLI arguments handling
	{
		CLI::App app{ "Multi - KMC (MKMC) ver. " MKMC_VER };

		configureArguments(argc, argv, params, app);

		// Add -- separator before positionals
		try {
			(app).parse(argc, argv);
		}
		catch (const CLI::ParseError& e) {
			std::string helpText;
			if (e.get_name() == "CallForHelp")
				helpText = app.help();
			else if (e.get_name() == "CallForAllHelp")
				helpText = app.help("", CLI::AppFormatMode::All);
			else
				return (app).exit(e);

			// Replace string
			size_t startPos = helpText.find("[OPTIONS]");
			if (startPos != std::string::npos)
				helpText.replace(startPos, std::string("[OPTIONS]").length(), "[OPTIONS] --");
			std::cout << helpText;
			return e.get_exit_code();
		}

		const bool CLIErrors = checkAndPrintArgumentsErrors(params);
		if (CLIErrors)
			std::exit(1);
	}

	if (params.mkmcParams.verbosity_level > 0)
		Logger::Inst().Enable();

	if (!params.readAdditionalParamsFromFiles())
		std::exit(1);

	const bool CLIAndParamsFromFilesErrors = checkAndPrintParamsFromFilesErrors(params);
	if (CLIAndParamsFromFilesErrors)
		std::exit(1);

	bool warningPrinted = checkAndPrintArgumentsWarnings(params);

	params.generateTempAndOutputFilesNames();
	warningPrinted |= params.adjustKMCPerformanceParams();
	warningPrinted |= params.adjustAnotherParams();
	params.readPhenotypes();

	Start start(params);
	Finish finish(params);

	try
	{
		Timer sequence_filter_init, kmc_timer, merger_timer, statistics_timer;

		bool wp = false;
		if (!start.verifyFiles(wp))
		{
			finish.finishProcessing();
			std::exit(1);
		}
		warningPrinted |= wp;

		if (warningPrinted) // any warning printed; insert distance before start stages
			std::cerr << std::endl;

		if (params.filterParams.filterKmersSequences)
		{
			SequenceFilterInit kmersFilter(params, sequence_filter_init);
			kmersFilter.prepareKmersSequencesToFilter();
		}

		std::cerr << "Starting k-mer counting..." << std::endl << std::endl;
		KMCRunner kmcRunner(params);
		kmc_timer.startTimer();
		kmcRunner.runKMCParallel();
		kmc_timer.stopTimer();

		MergerRunner dump_runner(params, merger_timer);
		DispatchKmerSize(params.stage1Params.GetKmerLen(), dump_runner);

		if (params.statisticsParams.generateNormalization || params.statisticsParams.generateEntropy || !params.statisticsParams.classificationMethods.empty())
		{
			std::vector<std::string> tasks;
			if (params.statisticsParams.generateNormalization)
				tasks.push_back("normalizing");
			if (!params.statisticsParams.correlationMethods.empty())
				tasks.push_back("computing correlation");
			if (params.statisticsParams.cvParams.cv)
				tasks.push_back("performing cross-validation");
			if (!params.statisticsParams.classificationMethods.empty())
				tasks.push_back("performing differential k-mers analysis");
			if (params.statisticsParams.generateEntropy)
				tasks.push_back("generating entropy");
			if (params.statisticsParams.runUMAP || params.statisticsParams.runPCA)
				tasks.push_back("reducing number of dimensions");

			std::cerr << "Starting " << MessagesUtilities::generateStartingSentence(tasks) << "..." << std::endl << std::endl;

			StatisticsGenerator statisticsGenerator(params);
			statistics_timer.startTimer();
			statisticsGenerator.generateStatisticsParallel();
			statistics_timer.stopTimer();
		}

		finish.finishProcessing();

		if (params.mutableParams.createdFastaFile)
		{
			std::cerr << "Preparing temporary FASTA file for sequences filtering out:\n";
			std::cerr << "\tStart: " << sequence_filter_init.getStartTime() << "\n";
			std::cerr << "\tEnd:   " << sequence_filter_init.getStopTime() << "\n";
		}
		std::cerr << "k-mer counting:\n";
		std::cerr << "\tStart: " << kmc_timer.getStartTime() << "\n";
		std::cerr << "\tEnd:   " << kmc_timer.getStopTime() << "\n";
		std::cerr << "Merging and dumping:\n";
		std::cerr << "\tStart: " << merger_timer.getStartTime() << "\n";
		std::cerr << "\tEnd:   " << merger_timer.getStopTime() << "\n";

		if (params.statisticsParams.generateNormalization)
		{
			if (params.statisticsParams.normalizationLearningWasSupplemented) // true in a case StatisticsGenerator detected, than DESeq2 learning data is missing
				std::cerr << "DESeq2 learning, normalizing and computing correlation:\n";
			else
				std::cerr << "Normalizing and computing correlation:\n";
			std::cerr << "\tStart: " << statistics_timer.getStartTime() << "\n";
			std::cerr << "\tEnd:   " << statistics_timer.getStopTime() << "\n";
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << e.what() << std::endl;
		finish.finishProcessing();
		std::exit(1);
	}
}
