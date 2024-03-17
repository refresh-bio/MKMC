#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include "kmc_core/kmc_runner.h"
#include "kmc_api/kmc_file.h"
#include "kmc_api/kmer_api.h"
#include "parameters.h"
#include "KMCRunner.h"
#include "KMCToolsRunner.h"
#include "Dump.h"
#include "Finish.h"
#include "Start.h"
#include "Version.h"
#include "time.hpp"
#include "FileGenerators.h"
#include "TasksFiller.h"
#include "Logger.h"


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
		<< "    E.g. -thr_rat0.5 and -thr2 mean that k-mers appearing at least twice in at least a half of the input files will be dumped.\n\n"
		<< "As `[options]` you can also pass optional parameters:\n"
		<< "  -k<len> - k-mer length (<len> from " << KMC::CfgConsts::min_k << " to " << KMC::CfgConsts::max_k << "; default: 25)\n"
		<< "  -f<a/q/m> - input in FASTA format (-fa), FASTQ format (-fq), or multi FASTA (-fm); mixing files is not supported (default: FASTQ)\n"
		<< "  -of<a,matrix> - output in FASTA format (-ofa) or matrix (-ofmatrix) (default: matrix)\n"
		<< "  -b - turn off transformation of k-mers into canonical form\n"
		<< "  -ci<X> - exclude k-mers occurring less than <X> times (if k-mer occurs less than <X> times in a file, it gets counter 0, but for this file only) (default: 1)\n"
		<< "  -cx<X> - exclude counting k-mers occurring more of than <X> times (if k-mer occurs more than <X> times in a file, it gets counter 0, but for this file only) (default: 4e9)\n"
		<< "  -cs<X> - maximal value of a counter (default: 65535)\n"
		<< "  -t<value> - number of threads (default: no. of CPU cores)\n"
		<< "  -wrk<X> - number of parallel k-mer counting tasks (default: 8)\n"
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

	for (uint32_t tmp_database_id = 0; tmp_database_id < mkmcParams.inputFilesPerSample.size(); ++tmp_database_id)
	{
		std::ostringstream sstreamKMCDir, sstreamKMC, sstreamTools;
		sstreamKMCDir << mkmcParams.tmpPath;
		sstreamKMC << mkmcParams.tmpPath;
		sstreamTools << mkmcParams.tmpPath;
		if (mkmcParams.tmpPath.back() != '/' && mkmcParams.tmpPath.back() != '\\')
		{
			sstreamKMCDir << "/";
			sstreamKMC << "/";
			sstreamTools << "/";
		}
		sstreamKMCDir << "kmc_tmp_" << std::setfill('0') << std::setw(5) << tmp_database_id;
		sstreamKMC << "kmc_db_" << std::setfill('0') << std::setw(5) << tmp_database_id;
		sstreamTools << "tools_db_" << std::setfill('0') << std::setw(5) << tmp_database_id;

		mkmcParams.kmcTmpDirs.push_back(sstreamKMCDir.str());
		mkmcParams.kmcOutputFiles.push_back(sstreamKMC.str());
		mkmcParams.toolsOutputFiles.push_back(sstreamTools.str());
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
	int i;

	bool was_m = false;
	bool was_r = false;

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
				std::cerr << "Error: Filtering threshold -thr should be a natural number.\n" << std::endl;
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
			if (*strEnd != '\0' || nThreads < 0)
			{
				std::cerr << "Error: Number of threads -t should be a natural number.\n" << std::endl;
				return false;
			}
			mkmcParams.nThreads = static_cast<uint32_t>(nThreads);
		}
		// Number of parallel KMC runs
		else if (strncmp(argv[i], "-wrk", 4) == 0)
		{
			int32_t nKMCWorkers = std::strtol(&argv[i][4], &strEnd, 10);
			if (*strEnd != '\0' || nKMCWorkers < 0)
			{
				std::cerr << "Error: Number of parallel k-mer counting tasks -wrk should be a natural number.\n" << std::endl;
				return false;
			}
			mkmcParams.nKMCWorkers = static_cast<uint32_t>(nKMCWorkers);
			mkmcParams.nKMCWorkersUserSet = true;
		}
		// k-mer length
		else if (strncmp(argv[i], "-k", 2) == 0)
		{
			int32_t k = std::strtol(&argv[i][2], &strEnd, 10);
			if (*strEnd != '\0' || k < 0)
			{
				std::cerr << "Error: k-mer length should be a natural number.\n" << std::endl;
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
				std::cerr << "Error: Maximal amount of RAM should be a natural number and be at least 2.\n" << std::endl;
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
				std::cerr << "Error: Minimal threshold of k-mers counter should be a natural number and be at least 1.\n" << std::endl;
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
				std::cerr << "Error: Maximal threshold of k-mers counter should be a natural number and be at least 1.\n" << std::endl;
				return false;
			}
			stage2Params.SetCutoffMax(static_cast<uint64_t>(cx));
		}
		// Maximal counter value
		else if (strncmp(argv[i], "-cs", 3) == 0)
		{
			int32_t cs = std::strtol(&argv[i][3], &strEnd, 10);
			if (*strEnd != '\0' || cs < 2)
			{
				std::cerr << "Error: Maximal k-mers counter value should be a natural number and be at least 2.\n" << std::endl;
				return false;
			}
			stage2Params.SetCounterMax(static_cast<uint64_t>(cs));
		}
		//output type
		else if (strncmp(argv[i], "-o", 2) == 0)
		{
			if (strncmp(argv[i] + 2, "fa", 2) == 0)
				mkmcParams.outputFileType = OutputFileType::FASTA;
			else if (strncmp(argv[i] + 2, "matrix", 6) == 0)
				mkmcParams.outputFileType = OutputFileType::Matrix;
			else
			{
				std::cerr << "Error: unsupported output type: " << argv[i] << " (use -ofa or -omatrix)\n" << std::endl;
				return false;
			}
		}
		// input type
		else if (strncmp(argv[i], "-f", 2) == 0)
		{
			if (strncmp(argv[i] + 2, "a", 1) == 0)
				stage1Params.SetInputFileType(KMC::InputFileType::FASTA);
			else if (strncmp(argv[i] + 2, "q", 1) == 0)
				stage1Params.SetInputFileType(KMC::InputFileType::FASTQ);
			else if (strncmp(argv[i] + 2, "m", 1) == 0)
				stage1Params.SetInputFileType(KMC::InputFileType::MULTILINE_FASTA);
			else
			{
				std::cerr << "Error: unsupported input type: " << argv[i] << " (use -fa, -fq, or -fm).\n" << std::endl;
				return false;
			}
		}
		// split dump file
		else if (strncmp(argv[i], "-dmp", 4) == 0)
		{
			mkmcParams.dumpStepSize = atoll(&argv[i][4]);
			mkmcParams.splitDumpOutput = true;
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

	std::string input_file_name = std::string(argv[i++]);

	mkmcParams.outputFile = argv[i++];

	mkmcParams.tmpPath = argv[i++];

	std::vector<std::string> input_file_names;
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

	fill_temporary_kmc_databases_names(params);

	return true;
}

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

		Timer kmc_timer, tools_timer, dump_timer;

		Start start(params);
		start.verifyFiles();

		std::cerr << "Starting k-mer counting..." << std::endl;
		KMCRunner kmcRunner(params);
		kmc_timer.startTimer();
		kmcRunner.runKMCParallel();
		kmc_timer.stopTimer();

		std::cerr << "\nStarting k-mer databases converting..." << std::endl;
		KMCToolsRunner kmcToolsRunner(params);
		tools_timer.startTimer();
		kmcToolsRunner.runKMCToolsParallel();
		tools_timer.stopTimer();

		Dump dump(params);
		std::cerr << "\nStarting dumping to file " << params.mkmcParams.outputFile << "..." << std::endl;
		if (params.mkmcParams.outputFileType == OutputFileType::Matrix)
		{
			dump_timer.startTimer();
			dump.dumpToFileParallel<MatrixFileGenerator>();
			dump_timer.stopTimer();
		}
		else if (params.mkmcParams.outputFileType == OutputFileType::FASTA)
		{
			dump_timer.startTimer();
			dump.dumpToFileParallel<FASTAFileGenerator>();
			dump_timer.stopTimer();
		}

		Finish finish(params);
		finish.finishProcessing();

		std::cerr << "\nKMC: \n";
		std::cerr << "\tStart: " << kmc_timer.getStartTime() << "\n";
		std::cerr << "\tEnd:   " << kmc_timer.getStopTime() << "\n";
		std::cerr << "KMC tools: \n";
		std::cerr << "\tStart: " << tools_timer.getStartTime() << "\n";
		std::cerr << "\tEnd:   " << tools_timer.getStopTime() << "\n";
		std::cerr << "Dump: \n";
		std::cerr << "\tStart: " << dump_timer.getStartTime() << "\n";
		std::cerr << "\tEnd:   " << dump_timer.getStopTime() << "\n";
	}
	catch (const std::exception& e)
	{
		std::cerr << e.what() << std::endl;
	}
}