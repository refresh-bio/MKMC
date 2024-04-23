#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <filesystem>
#include <limits>
#include "kmc_core/kmc_runner.h"
#include "kmc_api/kmc_file.h"
#include "kmc_api/kmer_api.h"
#include "parameters.h"
#include "KMCRunner.h"
#include "Dump.h"
#include "Finish.h"
#include "Start.h"
#include "Version.h"
#include "time.hpp"
#include "FileGenerators.h"
#include "TasksFiller.h"
#include "SequenceFilterInit.h"
#include "Logger.h"
#include "StatisticsGenerator.h"


//----------------------------------------------------------------------------------
// Show execution options of the software
void usage()
{
	std::cout << "Multi-KMC (MKMC) ver. " << MKMC_VER << "\n\n"
		<< "To run MKMC on Linux, type:\n"
		<< "  ./mkmc -k20 -thr_rat0.5 @<input_files> <output_file> <temp_dir>\n"
		<< "It will generate a matrix of 20-mers occurring in at least a half of the input files.\n\n"
		<< "General usage:\n"
		<< "  mkmc [options] @<input_files> <output_file> <temp_dir>\n\n"
		<< "The main parameters are as follows:\n"
		<< "  <input_files> - file with a list of samples names with input files names in specified (-f<a/q/m> switch) format (gzipped or not)\n"
		<< "  File example:\n"
		<< "    killifishretina1 kfA_1.fastq.gz kfA_2.fastq.gz\n"
		<< "    killifishretina2 kfB.fastq.gz\n"
		<< "  <output_file> - file where the matrix of k-mers counts or FASTA file will be dumped\n"
		<< "  <temp_dir> - a directory where temporary files will be stored\n"
		<< "  -thr<X> - filter out k-mers occuring less than <X> times... (default: 1)\n"
		<< "  -thr_rat<Y> - ... in a ratio <Y> of the input files (per k-mer sequence filtering) (default: 0.0)\n"
		<< "  -flt<X> - keep k-mers present in <X> file (FASTA or a set of the k-mers, one in each line) only\n"
		<< "  -n<freq/q> - generate normalized counts (frequency count/quantile normalization)\n"
		<< "  -p<X> - set a phenotype file <X> (a set of the integers, one in each line)\n"
		<< "  -cor<kendall/pearson/spearman> - determine correlations basing on a phenotype file (requires also -n and -p) (Kendall Tau/Pearson/Spearman correlation)\n"
		<< "    E.g. -thr_rat0.5 and -thr2 mean that k-mers appearing at least twice in at least a half of the input files will be dumped.\n\n"
		<< "As `[options]` you can also pass optional parameters:\n"
		<< "  -k<len> - k-mer length (<len> from " << KMC::CfgConsts::min_k << " to " << KMC::CfgConsts::max_k << "; default: 25)\n"
		<< "  -f<a/q/m> - input in FASTA format (-fa), FASTQ format (-fq), or multi FASTA (-fm); mixing files is not supported (default: FASTQ)\n"
		<< "  -of<a,matrix> - output in FASTA format (-ofa) or matrix (-ofmatrix) (default: matrix)\n"
		<< "  -on<X> - number of output files, reduce carefully (default: 512)\n"
		<< "  -b - turn off transformation of k-mers into canonical form\n"
		<< "  -ci<X> - exclude k-mers occurring less than <X> times (if k-mer occurs less than <X> times in a file, it gets counter 0, but for this file only) (default: 1)\n"
		<< "  -cx<X> - exclude counting k-mers occurring more of than <X> times (if k-mer occurs more than <X> times in a file, it gets counter 0, but for this file only) (default: 4e9)\n"
		<< "  -cs<X> - maximal value of a counter (default: 65535)\n"
		<< "  -t<value> - number of threads (default: 16 or no. of CPU cores)\n"
		<< "  -wrk<X> - number of parallel k-mer counting tasks (default: 4)\n"
		<< "  -m<X> - max amount of RAM in GB (from 2 to 1024); practically works only if -r is not set (default: 16)\n"
		<< "  -r - RAM only mode for k-mer counting\n"
		<< "  -v - verbose mode, shows progress\n"
		<< std::endl;
}

//----------------------------------------------------------------------------------
// Check if --help or --version was used
bool help_or_version(int argc, char** argv)
{
	const std::string version = "--version";
	const std::string help = "--help";
	for (int i = 1; i < argc; ++i)
	{
		if (argv[i] == version || argv[i] == help)
			return true;
	}
	return false;
}

void fill_temporary_kmc_databases_names(Params& params)
{
	MKMCParams& mkmcParams = params.mkmcParams;

	std::string tmpFilesTemplate = mkmcParams.tmpPath;
	if (mkmcParams.tmpPath.back() != '/' && mkmcParams.tmpPath.back() != '\\')
	{
		tmpFilesTemplate += static_cast<char>(std::filesystem::path::preferred_separator);
	}

	for (uint32_t tmp_database_id = 0; tmp_database_id < mkmcParams.inputFilesPerSample.size(); ++tmp_database_id)
	{
		std::ostringstream sstreamKMCDir, sstreamKMC;
		sstreamKMCDir << tmpFilesTemplate << "kmc_tmp_" << std::setfill('0') << std::setw(5) << tmp_database_id;
		sstreamKMC << tmpFilesTemplate << "kmc_db_" << std::setfill('0') << std::setw(5) << tmp_database_id;

		mkmcParams.kmcTmpDirs.push_back(sstreamKMCDir.str());
		mkmcParams.kmcOutputFiles.push_back(sstreamKMC.str());
	}

	params.statisticsParams.normFrequencyFileTmp = tmpFilesTemplate + params.statisticsParams.normFrequencyFileTmp;
	params.statisticsParams.normQuantileFileTmp = tmpFilesTemplate + params.statisticsParams.normQuantileFileTmp;

	const uint32_t nBinsDigits = static_cast<uint32_t>(std::log10(static_cast<double>(params.stage1Params.GetNBins()))) + 1;
	for (uint32_t binId = 0; binId < params.stage1Params.GetNBins(); ++binId) {
		std::ostringstream sstreamOutput;
		sstreamOutput << std::setfill('0') << std::setw(nBinsDigits) << binId;
		const std::string binIdStr = sstreamOutput.str();

		mkmcParams.outputMatrixFiles.push_back(mkmcParams.outputFilesTemplate + "_matrix_" + binIdStr);
		mkmcParams.outputFASTAFiles.push_back(mkmcParams.outputFilesTemplate + + "_" + binIdStr + ".fa");

		mkmcParams.outputFilesNormFrequency.push_back(mkmcParams.outputFilesTemplate + "_norm_" + binIdStr);
		mkmcParams.outputFilesNormQuantile.push_back(mkmcParams.outputFilesTemplate + "_norm_" + binIdStr);

		mkmcParams.outputFilesPearson.push_back(mkmcParams.outputFilesTemplate + "_pearson_" + binIdStr);
		mkmcParams.outputFilesSpearman.push_back(mkmcParams.outputFilesTemplate + "_spearman_" + binIdStr);
		mkmcParams.outputFilesKendall.push_back(mkmcParams.outputFilesTemplate + "_kendall_tau_" + binIdStr);
	}
}

//----------------------------------------------------------------------------------
// Parse the parameters
bool parse_parameters(int argc, char* argv[], Params& params)
{
	KMC::Stage1Params& stage1Params = params.stage1Params;
	KMC::Stage2Params& stage2Params = params.stage2Params;
	MKMCParams& mkmcParams = params.mkmcParams;
	FilterParams& filterParams = params.filterParams;
	StatisticsParams& statisticsParams = params.statisticsParams;
	int i;

	bool was_m = false;
	bool was_r = false;

	std::vector<OutputFileType> outputFileTypes;

	if (argc < 4)
		return false;

	for (i = 1; i < argc; ++i)
	{
		char* strEnd = nullptr;
		if (argv[i][0] != '-')
			break;
		// Filtering ratio
		if (strncmp(argv[i], "-thr_rat", 8) == 0) // must be before -t
		{
			filterParams.minKmersAboveThresholdRatio = std::strtof(&argv[i][8], &strEnd);
			if (*strEnd != '\0' || filterParams.minKmersAboveThresholdRatio < 0.0 || filterParams.minKmersAboveThresholdRatio > 1.0)
			{
				std::cerr << "Error: Filtering threshold -thr_rat should be a real number from a range [0, 1].\n" << std::endl;
				return false;
			}
		}
		// Filtering threshold
		else if (strncmp(argv[i], "-thr", 4) == 0)
		{
			int32_t minCountThreshold = std::strtol(&argv[i][4], &strEnd, 10);
			if (*strEnd != '\0' || minCountThreshold < 0)
			{
				std::cerr << "Error: Filtering threshold -thr should be a non-negative integer.\n" << std::endl;
				return false;
			}
			filterParams.minCountThreshold = static_cast<uint32_t>(minCountThreshold);
		}
		else if (strcmp(argv[i], "-keep") == 0)
		{
			mkmcParams.keepTmpFiles = true;
		}
		// Number of threads
		else if (strncmp(argv[i], "-t", 2) == 0)
		{
			int32_t nThreads = std::strtol(&argv[i][2], &strEnd, 10);
			if (*strEnd != '\0' || nThreads < 1)
			{
				std::cerr << "Error: Number of threads -t should be a positive integer.\n" << std::endl;
				return false;
			}
			mkmcParams.nThreads = static_cast<uint32_t>(nThreads);
		}
		// Number of parallel KMC runs
		else if (strncmp(argv[i], "-wrk", 4) == 0)
		{
			int32_t nKMCWorkers = std::strtol(&argv[i][4], &strEnd, 10);
			if (*strEnd != '\0' || nKMCWorkers < 1)
			{
				std::cerr << "Error: Number of parallel k-mer counting tasks -wrk should be a positive integer.\n" << std::endl;
				return false;
			}
			mkmcParams.nKMCWorkers = static_cast<uint32_t>(nKMCWorkers);
			mkmcParams.nKMCWorkersUserSet = true;
		}
		// k-mer length
		else if (strncmp(argv[i], "-k", 2) == 0)
		{
			int32_t k = std::strtol(&argv[i][2], &strEnd, 10);
			if (*strEnd != '\0' || k < 1)
			{
				std::cerr << "Error: k-mer length should be a positive integer.\n" << std::endl;
				return false;
			}
			stage1Params.SetKmerLen(static_cast<uint32_t>(k));
		}
		// Memory limit
		else if (strncmp(argv[i], "-m", 2) == 0)
		{
			int32_t maxRamGB = std::strtol(&argv[i][2], &strEnd, 10);
			if (*strEnd != '\0' || maxRamGB < 2)
			{
				std::cerr << "Error: Maximal amount of RAM should be a positive integer number and be at least 2.\n" << std::endl;
				return false;
			}
			mkmcParams.maxRamGB = static_cast<uint32_t>(maxRamGB);
			was_m = true;
		}
		// Minimum counter threshold
		else if (strncmp(argv[i], "-ci", 3) == 0)
		{
			int32_t ci = std::strtol(&argv[i][3], &strEnd, 10);
			if (*strEnd != '\0' || ci < 1)
			{
				std::cerr << "Error: Minimal threshold of k-mers counter should be a positive integer.\n" << std::endl;
				return false;
			}
			stage2Params.SetCutoffMin(static_cast<uint64_t>(ci));
		}
		// Maximum counter threshold
		else if (strncmp(argv[i], "-cx", 3) == 0)
		{
			int32_t cx = std::strtol(&argv[i][3], &strEnd, 10);
			if (*strEnd != '\0' || cx < 1)
			{
				std::cerr << "Error: Maximal threshold of k-mers counter should be a positive integer.\n" << std::endl;
				return false;
			}
			stage2Params.SetCutoffMax(static_cast<uint64_t>(cx));
		}
		// Maximal counter value
		else if (strncmp(argv[i], "-cs", 3) == 0)
		{
			int64_t cs = std::strtoll(&argv[i][3], &strEnd, 10);
			if (*strEnd != '\0' || cs < 2 || cs > static_cast<int64_t>(std::numeric_limits<uint32_t>::max()))
			{
				std::cerr << "Error: Maximal k-mers counter value should be a positivie integer and be at least 2.\n" << std::endl;
				return false;
			}
			stage2Params.SetCounterMax(static_cast<uint64_t>(cs));
		}
		// Number of KMC bins and output files
		else if (strncmp(argv[i], "-on", 3) == 0) // must be before -o
		{
			int32_t nKMCBins = std::strtol(&argv[i][3], &strEnd, 10);
			if (*strEnd != '\0' || nKMCBins < 1)
			{
				std::cerr << "Error: Number of output files should be a positive integer.\n" << std::endl;
				return false;
			}
			mkmcParams.nKMCBins = static_cast<uint32_t>(nKMCBins);
		}
		//output type
		else if (strncmp(argv[i], "-o", 2) == 0)
		{
			bool wasDuplication = false;
			OutputFileType outputFileType = OutputFileType::Matrix;

			if (strncmp(argv[i] + 2, "fa", 2) == 0)
				outputFileType = OutputFileType::FASTA;
			else if (strncmp(argv[i] + 2, "matrix", 6) == 0)
				outputFileType = OutputFileType::Matrix;
			else
			{
				std::cerr << "Error: unsupported output type: " << argv[i] << " (use -ofa or -omatrix)\n" << std::endl;
				return false;
			}

			for (auto fileType : outputFileTypes)
			{
				if (fileType == outputFileType)
				{
					std::cerr << "Warning: output format flag " << argv[i] << " was given multiple times." << std::endl;
					wasDuplication = true;
				}
			}
			if (!wasDuplication)
				outputFileTypes.push_back(outputFileType);
		}
		// File with k-mers to be filtered out
		else if (strncmp(argv[i], "-flt", 4) == 0) //  must be before -f
		{
			std::string fileName = &argv[i][4];
			if (fileName.empty())
			{
				std::cerr << "Error: file name of k-mers to be filtered out is not specified.\n" << std::endl;
				return false;
			}
			filterParams.filterKmersSequences = true;
			filterParams.inputKmersSequencesToFilterOut = fileName;
		}
		// input type
		else if (strncmp(argv[i], "-f", 2) == 0)
		{
			if (strncmp(argv[i] + 2, "a", 1) == 0)
				mkmcParams.inputFileType = KMC::InputFileType::FASTA;
			else if (strncmp(argv[i] + 2, "q", 1) == 0)
				mkmcParams.inputFileType = KMC::InputFileType::FASTQ;
			else if (strncmp(argv[i] + 2, "m", 1) == 0)
				mkmcParams.inputFileType = KMC::InputFileType::MULTILINE_FASTA;
			else
			{
				std::cerr << "Error: unsupported input type: " << argv[i] << " (use -fa, -fq, or -fm).\n" << std::endl;
				return false;
			}
		}
		// Generate statistics
		else if (strncmp(argv[i], "-n", 2) == 0)
		{
			std::string typeName = &argv[i][2];
			if (typeName == "freq")
			{
				statisticsParams.generateNormalization = true;
				statisticsParams.normalizationMethod = StatisticsParams::NormalizationMethod::frequency_count;
			}
			else if (typeName == "q")
			{
				statisticsParams.generateNormalization = true;
				statisticsParams.normalizationMethod = StatisticsParams::NormalizationMethod::quantile;
			}
			else
			{
				std::cerr << "Error: unsupported normalization method: " << typeName << " (use -nfreq, or -nq).\n" << std::endl;
				return false;
			}
		}
		else if (strncmp(argv[i], "-p", 2) == 0)
		{
			std::string fileName = &argv[i][2];
			if (fileName.empty())
			{
				std::cerr << "Error: phenotype file name is not specified.\n" << std::endl;
				return false;
			}
			statisticsParams.phenotypeFile = fileName;
		}
		else if (strncmp(argv[i], "-cor", 4) == 0)
		{
			std::string typeName = &argv[i][4];
			bool wasDuplication = false;
			StatisticsParams::CorrelationMethod correlationMethod;
			if (typeName == "kendall")
			{
				correlationMethod = StatisticsParams::CorrelationMethod::Kendall;
			}
			else if (typeName == "pearson")
			{
				correlationMethod = StatisticsParams::CorrelationMethod::Pearson;
			}
			else if (typeName == "spearman")
			{
				correlationMethod = StatisticsParams::CorrelationMethod::Spearman;
			}
			else
			{
				std::cerr << "Error: unsupported correlation method: " << typeName << " (use -corkendall, -corpearson, or -corspearman).\n" << std::endl;
				return false;
			}

			for (auto method : params.statisticsParams.correlationMethods)
			{
				if (method == correlationMethod)
				{
					std::cerr << "Warning: correlation method " << typeName << " was given multiple times." << std::endl;
					wasDuplication = true;
				}
			}
			if (!wasDuplication)
				params.statisticsParams.correlationMethods.push_back(correlationMethod);
		}
		else if (strncmp(argv[i], "-v", 2) == 0)
		{
			mkmcParams.verbosity_level++;
		}
		else if (strncmp(argv[i], "-r", 2) == 0)
		{
			stage1Params.SetRamOnlyMode(true);
			was_r = true;
		}
		else if (strncmp(argv[i], "-b", 2) == 0)
			stage1Params.SetCanonicalKmers(false);
	}

	if (argc - i < 3)
		return false;

	if (was_m && was_r)
	{
		std::cerr << "Warning: when -r parameter is given, limit specified with -m may be exceeded." << std::endl;
	}

	if (!outputFileTypes.empty())
		params.mkmcParams.outputFileTypes = outputFileTypes;

	if (!statisticsParams.correlationMethods.empty() && (!statisticsParams.generateNormalization || statisticsParams.phenotypeFile.empty()))
	{
		std::cerr << "Error: -cor parameter requires also -n and -p parameters.\n" << std::endl;
		return false;
	}

	if (statisticsParams.correlationMethods.empty() && !statisticsParams.phenotypeFile.empty())
	{
		std::cerr << "Warning: phenotype file was given (-p parameter), but no correlation -cor method is chosen; it wille be ignored." << std::endl;
	}

	std::string input_file_name = std::string(argv[i++]);

	mkmcParams.outputFilesTemplate = argv[i++];

	mkmcParams.tmpPath = argv[i++];

	if (input_file_name[0] != '@')
	{
		return false;
	}
	else
	{
		TasksFiller tasksFiller(mkmcParams, input_file_name);
		if (!tasksFiller.readSamples(params.mkmcParams.samples, params.mkmcParams.inputFilesPerSample))
			return false;
	}

	if (statisticsParams.generateNormalization)
	{
		bool generateMatrix = false;
		for (auto fileType : mkmcParams.outputFileTypes)
			if (fileType == OutputFileType::Matrix)
				generateMatrix = true;
		if (!generateMatrix)
		{
			std::cerr << "Warning: due to normalization generation, temporarily MKMC has to generate output matrix (-omatrix switch will be additionally applied)." << std::endl;
			mkmcParams.outputFileTypes.push_back(OutputFileType::Matrix);
		}
	}

	if (filterParams.filterKmersSequences) {
		std::ifstream in(filterParams.inputKmersSequencesToFilterOut);
		if (!in.good())
		{
			std::cerr << "Error: No " << params.filterParams.inputKmersSequencesToFilterOut << " file.\n" << std::endl;
			return false;
		}

		params.filterParams.kmersSequencesToFilterOutDB = params.mkmcParams.tmpPath + static_cast<char>(std::filesystem::path::preferred_separator) + "filter";
		params.mutableParams.kmersSequencesToFilterOut = params.mkmcParams.tmpPath + static_cast<char>(std::filesystem::path::preferred_separator) + "filter.fa";
	}

	if (!statisticsParams.correlationMethods.empty() && !statisticsParams.phenotypeFile.empty())
	{
		std::ifstream in(statisticsParams.phenotypeFile);
		if (!in.good())
		{
			std::cerr << "Error: No " << statisticsParams.phenotypeFile << " file.\n" << std::endl;
			return false;
		}
	}

	fill_temporary_kmc_databases_names(params);

	return true;
}

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
//----------------------------------------------------------------------------------
// Main function
int main(int argc, char** argv)
{
	if (argc == 1 || help_or_version(argc, argv))
	{
		usage();
		return 0;
	}

	try
	{
		Params params;
		if (!parse_parameters(argc, argv, params))
		{
			usage();
			return 0;
		}
		if (params.mkmcParams.verbosity_level > 0)
			Logger::Inst().Enable();

		params.setKMCParams();

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
			std::cerr << "\nStarting normalizing and correlation computing...\n";
			StatisticsGenerator statisticsGenerator(params);
			statistics_timer.startTimer();
			statisticsGenerator.generateStatisticsParallel();
			statistics_timer.stopTimer();
		}

		Finish finish(params);
		finish.finishProcessing();

		if (params.mutableParams.createdFastaFile)
		{
			std::cerr << "\nPreparing temporary FASTA file for sequences filtering out: \n";
			std::cerr << "\tStart: " << sequence_filter_init.getStartTime() << "\n";
			std::cerr << "\tEnd:   " << sequence_filter_init.getStopTime() << "\n";
		}
		std::cerr << "\nKMC: \n";
		std::cerr << "\tStart: " << kmc_timer.getStartTime() << "\n";
		std::cerr << "\tEnd:   " << kmc_timer.getStopTime() << "\n";
		std::cerr << "Dump: \n";
		std::cerr << "\tStart: " << dump_timer.getStartTime() << "\n";
		std::cerr << "\tEnd:   " << dump_timer.getStopTime() << "\n";

		if (params.statisticsParams.generateNormalization)
		{
			std::cerr << "Normalization and correlation: \n";
			std::cerr << "\tStart: " << statistics_timer.getStartTime() << "\n";
			std::cerr << "\tEnd:   " << statistics_timer.getStopTime() << "\n";
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << e.what() << std::endl;
	}
}