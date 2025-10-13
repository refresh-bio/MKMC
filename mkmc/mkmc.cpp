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

		Logger::Inst().Log(std::string("\nStarting merging samples and dumping to ") + (tasks.size() > 1 ? "files " : "file ") + MessagesUtilities::generateSentence(tasks) + "...");

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

	app.add_option("input_samples_file", mkmcParams.inputFileName, "file with a list of samples and input files in specified (-f parameter) format (gzipped or not)")->required()->check(CLI::ExistingFile);
	app.add_option("output_files_template", mkmcParams.outputFilesTemplate, "template (prefix) of output files names")->required();
	app.add_option("temp_dir", mkmcParams.tmpPath, "directory for temporary files")->required();

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
	filteringGroup->add_option_function("--flt", fltCallback, "keep k-mers present in a specified file (FASTA or a set of the k-mers, one per line) only; if -b is not set, the k-mers are converted to canonical form")->check(CLI::ExistingFile);

	CLI::Option_group* correlationGroup = app.add_option_group("correlation and normalization");

	std::map<std::string, StatisticsParams::NormalizationMethod> valuesMap{ {"deseq", StatisticsParams::NormalizationMethod::deseq2}, {"freq", StatisticsParams::NormalizationMethod::frequency_count }, {"q", StatisticsParams::NormalizationMethod::quantile } };
	std::function<void(const decltype(statisticsParams.normalizationMethod)&)> nCallback = [&](const decltype(statisticsParams.normalizationMethod)& normalizationMethod)
	{
		statisticsParams.normalizationMethod = normalizationMethod;
		statisticsParams.generateNormalization = true;
	};
	n = correlationGroup->add_option_function("-n", nCallback, "normalize counts (DESeq2/frequency count/quantile normalization) before use")->transform(CLI::CheckedTransformer(valuesMap, CLI::ignore_case));

	correlationGroup->add_flag("--save_n", statisticsParams.saveNormalization, "save normalized matrix to file")->needs(n);

	std::map<std::string, StatisticsParams::CorrelationMethod> correlationValuesMap{ {"pearson", StatisticsParams::CorrelationMethod::Pearson }, { "spearman", StatisticsParams::CorrelationMethod::Spearman }, {"kendall", StatisticsParams::CorrelationMethod::Kendall } };
	cor = correlationGroup->add_option("--cor", statisticsParams.correlationMethods, "compute correlation coefficients, basing on a phenotype file (Kendall Tau/Pearson/Spearman correlation)")->transform(CLI::CheckedTransformer(correlationValuesMap));

	std::function<void(const std::string&)> pCallback = [&](const std::string& fileName)
	{
		phenotypes.correlationPhenotype.setFileName(fileName);
	};
	p = correlationGroup->add_option_function("-p", pCallback, "set a phenotype file (a sequence of integers, one in each line)")->check(CLI::ExistingFile)->needs(cor);

	CLI::Option_group* diffGroup = app.add_option_group("differential k-mers analysis");
	typedef StatisticsParams::DifferentialAnalysisMethod DAMethod;
	std::map<std::string, DAMethod> differentialAnalysisValuesMap{ {"ttest", DAMethod::TTest }, {"snr", DAMethod::SNR }, {"wrs", DAMethod::WilcoxonRankSum }, {"dids", DAMethod::DIDS }, {"anova", DAMethod::ANOVA } };
	differentialAnalysis = diffGroup->add_option("--diff", statisticsParams.classificationMethods, "perform differential k-mers analysis (ANOVA, DIDS, Signal to Noise ratio, T-Test, Wilcoxon-rank sum (Mann-Whitney U test)); all except T-Test need -n; counts for T-Test are always unnormalized, increased by 1, and logarithmized")->transform(CLI::CheckedTransformer(differentialAnalysisValuesMap));

	typedef StatisticsParams::DifferentialAnalysisCorrectionMethod CorrectionMethod;
	std::map<std::string, StatisticsParams::DifferentialAnalysisCorrectionMethod> differentialAnalysisCorrectionValuesMap{ { "b", CorrectionMethod::Bonferroni }, { "hb", CorrectionMethod::HolmBonferroni }, { "bh", CorrectionMethod::BenjaminiHochberg }, { "by", CorrectionMethod::BenjaminiYekutieli } };
	std::function<void(const decltype(statisticsParams.classificationPValueCorrection)&)> pcorrCallback = [&](const decltype(statisticsParams.classificationPValueCorrection)& classificationPValueCorrection)
	{
		statisticsParams.classificationPValueCorrection = classificationPValueCorrection;
		statisticsParams.correctPvalues = true;
	};
	pvalCorr = diffGroup->add_option_function("--pval_corr", pcorrCallback, "correct p-values of differential k-mers analysis (Bonferroni, Benjamini-Hochberg, Benjamini-Yekutieli, Holm-Bonferroni); store statistically significant k-mers also in separated files; useful for ANOVA, T-Test, Wilcoxon-rank sum")->transform(CLI::CheckedTransformer(differentialAnalysisCorrectionValuesMap))->needs(differentialAnalysis);

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

	auto umap = dimReductionGroup->add_flag("--umap", statisticsParams.runUMAP, "reduce dimensionality of normalized matrix with UMAP")->needs(n);
	auto pca = dimReductionGroup->add_flag("--pca", statisticsParams.runPCA, "reduce dimensionality of normalized matrix with PCA")->needs(n);

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
	optionalGroup->add_flag_callback("-r", rCallback, "count k-mers in RAM only");

	optionalGroup->add_flag("-v", mkmcParams.verbosity_level, "verbose mode, shows progress and minor warnings");

	CLI::Option_group* debugGroup = app.add_option_group("debug parameters");
	auto keep = debugGroup->add_flag("--keep", mkmcParams.keepTmpFiles, "keep temporary files and binary results file");
	debugGroup->add_flag("--reuse-db", mkmcParams.reuseDBFiles, "reuse samples and filtering databases (if possible)");
	debugGroup->add_flag("--learn-deseq", statisticsParams.learnDeseq2, "collect data for DESeq2 normalization (not necessary for -n deseq, but useful for further --reuse-db)")->needs(keep);

	debugGroup->add_option("--on", mkmcParams.nKMCBins, "suggested number of internal bins, modify carefully")->check(CLI::PositiveNumber)->default_val(mkmcParams.nKMCBins);

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
		(statisticsParams.classificationMethods.size() > 1 ||
			(statisticsParams.classificationMethods.size() == 1 && statisticsParams.classificationMethods.front() != StatisticsParams::DifferentialAnalysisMethod::TTest)))
	{
		Logger::Inst().Log("Error: Differential analysis methods (except T-Test) require normalization (-n).");
		return true;
	}

	if (statisticsParams.nDimensionReductionUserDefined &&
		!(statisticsParams.runPCA || statisticsParams.runUMAP))
	{
		Logger::Inst().Log("Error: Number of dimensions (--dimensions) requires dimensionality reduction algorithm (--umap or --pca).");
		return true;
	}

	if (statisticsParams.correctPvalues &&
		!isDAMethod(StatisticsParams::DifferentialAnalysisMethod::ANOVA) &&
		!isDAMethod(StatisticsParams::DifferentialAnalysisMethod::TTest) &&
		!isDAMethod(StatisticsParams::DifferentialAnalysisMethod::WilcoxonRankSum))
	{
		Logger::Inst().Log("Error: --pval_corr requires differential k-mers analysis with ANOVA, T-Test, or Wilcoxon-rank sum (Mann-Whitney U test) (--diff).");
		return true;
	}

	if (statisticsParams.didsModeUserDefined && !isDAMethod(StatisticsParams::DifferentialAnalysisMethod::DIDS))
	{
		Logger::Inst().Log("Error: --dids-mode requires differential k-mers analysis with DIDS (--diff).");
		return true;
	}

	if (statisticsParams.nTopUserDefined &&
		statisticsParams.correlationMethods.empty() &&
		!statisticsParams.generateEntropy &&
		!isDAMethod(StatisticsParams::DifferentialAnalysisMethod::SNR) &&
		!isDAMethod(StatisticsParams::DifferentialAnalysisMethod::DIDS))
	{
		Logger::Inst().Log("Error: --n_top requires correlation (--cor) or entropy (--entropy) or differential k-mers analysis with SNR or DIDS (--diff).");
		return true;
	}

	if (statisticsParams.nTopUserDefined &&
		statisticsParams.nTop < 1)
	{
		Logger::Inst().Log("Error: --n_top has to be at least 1.");
		return true;
	}

	return false;
}



bool checkAndPrintDataFromFilesVsParamsErrors(const Params& params)
{
	if ((params.statisticsParams.runPCA || params.statisticsParams.runUMAP) &&
		(params.statisticsParams.nDimensionReduction < 1 || params.statisticsParams.nDimensionReduction >= params.mkmcParams.samples.size()))
	{
		Logger::Inst().Log("Error: Number of dimensions (--dimensions) must be at least 1 and lower than number of samples.");
		return true;
	}

	if (params.statisticsParams.cvParams.cv)
		if (params.statisticsParams.cvParams.p == 0 || params.statisticsParams.cvParams.p >= params.mkmcParams.samples.size() || params.mkmcParams.samples.size() % params.statisticsParams.cvParams.p != 0)
		{
			Logger::Inst().Log("Error: Number of samples to leave in cross-validation (--leave) has to be positive and be a factor of a number of samples.");
			return true;
		}
	return false;
}


bool checkAndPrintArgumentsWarnings(const Params& params)
{
	bool result = false;
	if (params.mkmcParams.samples.size() <= 8 && std::find(params.statisticsParams.classificationMethods.begin(), params.statisticsParams.classificationMethods.end(), StatisticsParams::DifferentialAnalysisMethod::WilcoxonRankSum) != params.statisticsParams.classificationMethods.end())
	{
		Logger::Inst().Log("Warning: Wilcoxon-rank sum (Mann-Whitney U test) uses approximate algorithm, thus for less than 9 samples its results may be slightly different than in e.g. SciPy.", 1);
		result = true;
	}

	if (params.statisticsParams.cvParams.seedUserDefined && params.statisticsParams.cvParams.p == 1)
	{
		Logger::Inst().Log("Warning: as --leave parameter is set to 1, LOOCV will be performed, which does not need randomness (--cv-seed parameter will be ignored).", 1);
		result = true;
	}

	if (params.statisticsParams.learnDeseq2 && params.statisticsParams.normalizationMethod == StatisticsParams::NormalizationMethod::deseq2)
	{
		Logger::Inst().Log("Warning: --learn-deseq is not necessary, when DESeq2 normalization is performed (-n deseq is set).", 1);
		result = true;
	}

	return result;
}



bool verifyDBsReusability(const Params& params)
{
	if (params.filterParams.filterKmersSequences)
	{
		try
		{
			kmcdb::MetadataReader sequencesToFilterMetadataReader(params.filterParams.kmersSequencesToFilterOutDB, false);
			if (params.stage1Params.GetKmerLen() != sequencesToFilterMetadataReader.GetConfig().kmer_len)
				return false;
		}
		catch (const std::runtime_error&)
		{
			return false;
		}
	}

	try
	{
		kmcdb::MetadataReader matrixMetadataReader(params.mkmcParams.outputMatrixBinFile, false);
		if (params.stage1Params.GetKmerLen() != matrixMetadataReader.GetConfig().kmer_len)
			return false;

		if (matrixMetadataReader.GetConfig().num_samples != params.mkmcParams.samples.size())
			return false;
	}
	catch (const std::runtime_error&)
	{
		return false;
	}

	if (params.statisticsParams.generateNormalization)
	{
		try
		{
			MatrixStatsReader statsReader(params.mkmcParams.normLearningBinFile);
			std::vector<uint8_t> tmp;
			// For DESeq2 there is possibility to supplement required learning data later, but here we verify reusability, thus data should be consistent (the more, --learn-deseq is available).
			// For another methods the learning data should always be present, if file exists.
			if (!statsReader.Get(StatisticsParams::getNormalizationMethodStreamName(params.statisticsParams.normalizationMethod), tmp))
				return false;
		}
		catch (const std::runtime_error&)
		{
			return false;
		}
	}

	// Actually, values of old --thr and --thr_rat parameters should be equal to current,
	// however currently it is impossible to compare them.

	return true;
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
			std::cout << helpText; // Print to cout, as it is CLI11 default behaviour.
			return e.get_exit_code();
		}

		const bool CLIErrors = checkAndPrintArgumentsErrors(params);
		if (CLIErrors)
			std::exit(1);
	}

	Logger::Inst().Enable(params.mkmcParams.verbosity_level, params.mkmcParams.getLogFileName());

	bool warningPrinted = false;
	bool filterMsgPrinted = false;

	if (!params.readAdditionalDataFromFiles(warningPrinted))
		std::exit(1);

	const bool CLIAndParamsFromFilesErrors = checkAndPrintDataFromFilesVsParamsErrors(params);
	if (CLIAndParamsFromFilesErrors)
		std::exit(1);

	warningPrinted |= checkAndPrintArgumentsWarnings(params);

	params.generateTempAndOutputFilesNames();
	warningPrinted |= params.adjustKMCPerformanceParams();
	warningPrinted |= params.adjustAnotherParams();

	const bool phenotypesReadSuccess = params.readPhenotypes();
	if (!phenotypesReadSuccess)
		std::exit(1);

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

		if (params.filterParams.filterKmersSequences)
		{
			if (warningPrinted) // any warning printed; insert distance before start stages
			{
				Logger::Inst().Log("", 1);
				warningPrinted = false;
			}

			SequenceFilterInit kmersFilter(params, sequence_filter_init);
			filterMsgPrinted = kmersFilter.prepareKmersSequencesToFilter();
		}

		bool dbsReusable = false;
		if (!params.mkmcParams.reuseDBFiles)
		{
			// any msg printed; insert distance before start stages
			if (filterMsgPrinted)
				Logger::Inst().Log("");
			else if (warningPrinted)
				Logger::Inst().Log("", 1);
			
			Logger::Inst().Log("Starting k-mer counting...");
		}
		else {
			dbsReusable = verifyDBsReusability(params);
			// any filter msg printed; insert distance after that stage
			if (filterMsgPrinted)
				Logger::Inst().Log("");
			if (!dbsReusable)
			{
				// any warning printed; insert distance before start stages
				if (!filterMsgPrinted && warningPrinted)
					Logger::Inst().Log("", 1);
				Logger::Inst().Log("Samples databases does not exist or are not possible to reuse. Starting k-mer counting...");
			}
			else
				Logger::Inst().Log("Info: Samples databases exist and outwardly seem to be possible to reuse.");
		}

		if (!dbsReusable)
		{
			KMCRunner kmcRunner(params);
			kmc_timer.startTimer();
			kmcRunner.runKMCParallel();
			kmc_timer.stopTimer();

			MergerRunner dump_runner(params, merger_timer);
			DispatchKmerSize(params.stage1Params.GetKmerLen(), dump_runner);
		}

		bool computeStatistics = false;
		if (params.statisticsParams.generateNormalization || params.statisticsParams.generateEntropy || !params.statisticsParams.classificationMethods.empty())
		{
			computeStatistics = true;

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

			Logger::Inst().Log("\nStarting " + MessagesUtilities::generateSentence(tasks) + "...");

			StatisticsGenerator statisticsGenerator(params);
			statistics_timer.startTimer();
			statisticsGenerator.generateStatisticsParallel();
			statistics_timer.stopTimer();
		}

		finish.finishProcessing();

		Logger::Inst().Log("");
		if (!dbsReusable)
		{
			if (params.mutableParams.createdFastaFile)
			{
				Logger::Inst().Log("Preparing k-mers for filtering:");
				Logger::Inst().Log("\tStart: " + sequence_filter_init.getStartTime());
				Logger::Inst().Log("\tEnd:   " + sequence_filter_init.getStopTime());
			}
			Logger::Inst().Log("k-mer counting:");
			Logger::Inst().Log("\tStart: " + kmc_timer.getStartTime());
			Logger::Inst().Log("\tEnd:   " + kmc_timer.getStopTime());
			Logger::Inst().Log("Merging and dumping:");
			Logger::Inst().Log("\tStart: " + merger_timer.getStartTime());
			Logger::Inst().Log("\tEnd:   " + merger_timer.getStopTime());
		}

		if (computeStatistics)
		{
			if (params.statisticsParams.normalizationLearningWasSupplemented) // true in a case StatisticsGenerator detected, than DESeq2 learning data is missing
			{
				std::vector<std::string> tasks{ "DESeq2 learning" };
				if (params.statisticsParams.generateNormalization)
					tasks.push_back("normalizing");
				if (!params.statisticsParams.correlationMethods.empty() || !params.statisticsParams.classificationMethods.empty() || params.statisticsParams.generateEntropy || params.statisticsParams.runUMAP || params.statisticsParams.runPCA)
					tasks.push_back("computing statistics");
				Logger::Inst().Log(MessagesUtilities::generateSentence(tasks) + ":");
			}
			else
			{
				std::vector<std::string> tasks;
				if (params.statisticsParams.generateNormalization)
					tasks.push_back("normalizing");
				if (!params.statisticsParams.correlationMethods.empty() || !params.statisticsParams.classificationMethods.empty() || params.statisticsParams.generateEntropy || params.statisticsParams.runUMAP || params.statisticsParams.runPCA)
					tasks.push_back("computing statistics");
				Logger::Inst().Log(MessagesUtilities::generateSentence(tasks, true) + ":");
			}
			Logger::Inst().Log("\tStart: " + statistics_timer.getStartTime());
			Logger::Inst().Log("\tEnd:   " + statistics_timer.getStopTime());
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << e.what() << std::endl;
		finish.finishProcessing();
		std::exit(1);
	}
}
