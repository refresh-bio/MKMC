#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <map>
#include <functional>
#include "kmc_core/kmc_runner.h"
#include "parameters.h"
#include "KMCRunner.h"
#include "Dump.h"
#include "Finish.h"
#include "Start.h"
#include "Version.h"
#include "time.hpp"
#include "SequenceFilterInit.h"
#include "Logger.h"
#include "StatisticsGenerator.h"
#include "CLI11/include/CLI/CLI.hpp"



class DumpRunner
{
	Params& params;
	Timer& dump_timer;
public:
	DumpRunner(Params& params, Timer& dump_timer):
		params(params),
		dump_timer(dump_timer)
	{

	}
	template<unsigned SIZE>
	void Run()
	{
		std::cerr << "\nStarting dumping to files " << params.mkmcParams.outputFilesTemplate << "_X..." << std::endl;
		Dump<SIZE> dump(params);
		dump_timer.startTimer();
		dump.dumpToFileParallel();
		dump_timer.stopTimer();
	}
};



void createArguments(int argc, char** argv, Params& params, CLI::App& app)
{
	DefaultKMCParams& defaultKMCParams = params.defaultKMCParams;
	KMC::Stage1Params& stage1Params = params.stage1Params;
	KMC::Stage2Params& stage2Params = params.stage2Params;
	MKMCParams& mkmcParams = params.mkmcParams;
	FilterParams& filterParams = params.filterParams;
	StatisticsParams& statisticsParams = params.statisticsParams;
	Phenotypes& phenotypes = params.phenotypes;

	// set KMC defaults
	stage1Params.SetKmerLen(defaultKMCParams.k);
	stage2Params.SetCutoffMin(defaultKMCParams.ci);
	stage2Params.SetCutoffMax(defaultKMCParams.cx);
	stage2Params.SetCounterMax(defaultKMCParams.cs);

	CLI::Option* p = nullptr, * n = nullptr, * cor = nullptr, * differentialAnalysis = nullptr, * c = nullptr;

	app.add_option("input_samples_file", mkmcParams.inputFileName, "file with a list of samples names with input files names in specified (-f parameter) format (gzipped or not)")->required()->check(CLI::ExistingFile);
	app.add_option("output_files", mkmcParams.outputFilesTemplate, "file where the matrix of k-mers counts or FASTA file will be dumped")->required();
	app.add_option("temp_dir", mkmcParams.tmpPath, "a directory where temporary files will be stored")->required()->check(CLI::ExistingDirectory);

	std::function<void(const uint32_t&)> kCallback = [&](const uint32_t& k)
	{
		stage1Params.SetKmerLen(k);
	};
	app.add_option_function("-k", kCallback, "k-mer length")->check(CLI::Range(KMC::CfgConsts::min_k, KMC::CfgConsts::max_k))->default_val(defaultKMCParams.k);

	app.add_option("--thr", filterParams.minCountThreshold, "filter out k-mers occuring less than specified number of times...")->check(CLI::PositiveNumber)->default_val(filterParams.minCountThreshold);
	app.add_option("--thr_rat", filterParams.minKmersAboveThresholdRatio, "... in a specified ratio of the input files (see example)")->check(CLI::Range(0.0, 1.0))->default_val(filterParams.minKmersAboveThresholdRatio);

	std::function<void(const decltype(filterParams.inputKmersSequencesToFilterOut)&)> fltCallback = [&](const decltype(filterParams.inputKmersSequencesToFilterOut)& fileName)
	{
		filterParams.inputKmersSequencesToFilterOut = fileName;
		filterParams.filterKmersSequences = true;
	};
	app.add_option_function("--flt", fltCallback, "keep k-mers present in a specified file (FASTA or a set of the k-mers, one in each line) only")->check(CLI::ExistingFile);

	std::map<std::string, StatisticsParams::NormalizationMethod> valuesMap{ {"freq", StatisticsParams::NormalizationMethod::frequency_count }, {"q", StatisticsParams::NormalizationMethod::quantile } };
	std::function<void(const decltype(statisticsParams.normalizationMethod)&)> nCallback = [&](const decltype(statisticsParams.normalizationMethod)& normalizationMethod)
	{
		statisticsParams.normalizationMethod = normalizationMethod;
		statisticsParams.generateNormalization = true;
	};
	n = app.add_option_function("-n", nCallback, "generate normalized counts (frequency count/quantile normalization)")->transform(CLI::CheckedTransformer(valuesMap, CLI::ignore_case));

	std::map<std::string, StatisticsParams::CorrelationMethod> correlationValuesMap{ {"pearson", StatisticsParams::CorrelationMethod::Pearson }, { "spearman", StatisticsParams::CorrelationMethod::Spearman }, {"kendall", StatisticsParams::CorrelationMethod::Kendall } };
	cor = app.add_option("--cor", statisticsParams.correlationMethods, "compute correlation cofficients with specified methods, basing on a phenotype file (Kendall Tau/Pearson/Spearman correlation)")->transform(CLI::CheckedTransformer(correlationValuesMap));

	std::function<void(const std::string&)> pCallback = [&](const std::string& fileName)
	{
		phenotypes.correlationPhenotype.setFileName(fileName);
	};
	p = app.add_option_function("-p", pCallback, "set a phenotype file (a set of the integers, one in each line)")->check(CLI::ExistingFile)->needs(cor);

	std::map<std::string, StatisticsParams::DifferentialAnalysisMethod> differentialAnalysisValuesMap{ {"t", StatisticsParams::DifferentialAnalysisMethod::TTest } };
	differentialAnalysis = app.add_option("--diff", statisticsParams.classificationMethods, "perform differential k-mers analysis (T-Test)")->transform(CLI::CheckedTransformer(differentialAnalysisValuesMap));

	std::function<void(const std::string&)> cCallback = [&](const std::string& fileName)
	{
		phenotypes.differentialAnalysisPhenotype.setFileName(fileName);
	};
	c = app.add_option_function("-c", cCallback, "set a phenotype file for differential k-mers analysis (a set of the natural numbers or text labels, one in each line)")->check(CLI::ExistingFile)->needs(differentialAnalysis);

	app.add_flag("--entropy", statisticsParams.generateEntropy, "generate k-mers counts entropy")->default_val(statisticsParams.generateEntropy);

	CLI::Option_group* optionalGroup = app.add_option_group("optional parameters");

	std::map<std::string, KMC::InputFileType> inputValuesMap{ {"fa", KMC::InputFileType::FASTA }, {"fq", KMC::InputFileType::FASTQ }, { "mf", KMC::InputFileType::MULTILINE_FASTA } };
	optionalGroup->add_option("-f", mkmcParams.inputFileType, "input format (FASTA, FASTQ or multi-FASTA); mixing files is not supported")->transform(CLI::CheckedTransformer(inputValuesMap, CLI::ignore_case))->default_val(mkmcParams.inputFileType)->default_str("fq");

	std::map<std::string, OutputFileType> outputValuesMap{ {"fa", OutputFileType::FASTA }, {"matrix", OutputFileType::Matrix } };
	optionalGroup->add_option("-o", mkmcParams.outputFileTypes, "output format (FASTA or matrix)")->transform(CLI::CheckedTransformer(outputValuesMap, CLI::ignore_case))->default_val(mkmcParams.outputFileTypes)->default_str("matrix");

	optionalGroup->add_option("--on", mkmcParams.nKMCBins, "number of output files, reduce carefully")->check(CLI::PositiveNumber)->default_val(mkmcParams.nKMCBins);

	std::function<void()> bCallback = [&]()
	{
		stage1Params.SetCanonicalKmers(false);
	};
	optionalGroup->add_flag_callback("-b", bCallback, "turn off transformation of k-mers into canonical form");

	std::function<void(const uint32_t&)> ciCallback = [&](const uint32_t& ci) // currently 32 bits
	{
		stage2Params.SetCutoffMin(static_cast<uint64_t>(ci));
	};
	optionalGroup->add_option_function("--ci", ciCallback, "exclude k-mers occurring less than specified number of times (if k-mer occurs less than --ci times in a sample, it gets counter 0, but for this sample only)")->check(CLI::PositiveNumber)->default_val(defaultKMCParams.ci);
	std::function<void(const uint32_t&)> cxCallback = [&](const uint32_t& cx) // currently 32 bits
	{
		stage2Params.SetCutoffMax(static_cast<uint64_t>(cx));
	};
	optionalGroup->add_option_function("--cx", cxCallback, "exclude counting k-mers occurring more than specified number of times (if k-mer occurs more than --cx times in a sample, it gets counter 0, but for this sample only)")->check(CLI::PositiveNumber)->default_val(static_cast<uint32_t>(defaultKMCParams.cx));

	std::function<void(const uint32_t&)> csCallback = [&](const uint32_t& cs) // currently 32 bits
	{
		stage2Params.SetCounterMax(static_cast<uint64_t>(cs));
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
	debugGroup->add_flag("--keep", mkmcParams.keepTmpFiles, "keep temporary files");

	cor->needs(n)->needs(p);
	differentialAnalysis->needs(c);

	app.footer("Example: to run MKMC, type:\n"
		"    ./mkmc -k 20 --thr_rat 0.5 input_files_list.txt output tmp\n"
		"It will generate a matrix of 20-mers occurring in at least a half of the input files.\n"
		"    ./mkmc -k 20 --thr 2 --thr_rat 0.5 input_files_list.txt output tmp\n"
		"It will generate a matrix of 20-mers occurring at least twice in at least a half of the input files.\n\n"
		"input_files_list.txt example:\n"
		"    killifishretina1 kfA_1.fastq.gz kfA_2.fastq.gz\n"
		"    killifishretina2 kfB.fastq.gz");
}



//----------------------------------------------------------------------------------
// Main function
int main(int argc, char** argv)
{
	Params params;
	CLI::App app{ "Multi - KMC (MKMC) ver. " MKMC_VER };

	createArguments(argc, argv, params, app);
	CLI11_PARSE(app, argc, argv);

	if (!params.readAdditionalParamsFromFiles())
	{
		std::exit(1);
	}

	params.generateTempAndOutputFilesNames();
	params.adjustKMCPerformanceParams();
	params.adjustAnotherParams();
	params.readPhenotypes();

	if (params.mkmcParams.verbosity_level > 0)
		Logger::Inst().Enable();

	try
	{
		Timer sequence_filter_init, kmc_timer, dump_timer, statistics_timer;

		Start start(params);
		start.verifyFiles();

		if (params.filterParams.filterKmersSequences)
		{
			SequenceFilterInit kmersFilter(params, sequence_filter_init);
			kmersFilter.prepareKmersSequencesToFilter();
		}

		if (params.mutableParams.createdFastaFile)
			std::cerr << '\n';
		std::cerr << "Starting k-mer counting..." << std::endl;
		KMCRunner kmcRunner(params);
		kmc_timer.startTimer();
		kmcRunner.runKMCParallel();
		kmc_timer.stopTimer();

		DumpRunner dump_runner(params, dump_timer);
		DispatchKmerSize(params.stage1Params.GetKmerLen(), dump_runner);

		if (params.statisticsParams.generateNormalization)
		{
			std::cerr << "\nStarting normalizing and computing correlation...\n";
			StatisticsGenerator statisticsGenerator(params);
			statistics_timer.startTimer();
			statisticsGenerator.generateStatisticsParallel();
			statistics_timer.stopTimer();
		}

		Finish finish(params);
		finish.finishProcessing();

		if (params.mutableParams.createdFastaFile)
		{
			std::cerr << "\nPreparing temporary FASTA file for sequences filtering out:\n";
			std::cerr << "\tStart: " << sequence_filter_init.getStartTime() << "\n";
			std::cerr << "\tEnd:   " << sequence_filter_init.getStopTime() << "\n";
		}
		std::cerr << "\nk-mer counting:\n";
		std::cerr << "\tStart: " << kmc_timer.getStartTime() << "\n";
		std::cerr << "\tEnd:   " << kmc_timer.getStopTime() << "\n";
		std::cerr << "Dumping:\n";
		std::cerr << "\tStart: " << dump_timer.getStartTime() << "\n";
		std::cerr << "\tEnd:   " << dump_timer.getStopTime() << "\n";

		if (params.statisticsParams.generateNormalization)
		{
			std::cerr << "Normalizing and computing correlation:\n";
			std::cerr << "\tStart: " << statistics_timer.getStartTime() << "\n";
			std::cerr << "\tEnd:   " << statistics_timer.getStopTime() << "\n";
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << e.what() << std::endl;
	}
}
