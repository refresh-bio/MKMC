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



//----------------------------------------------------------------------------------
// Show execution options of the software
void usage()
{
	std::cout << "K-Mer Counter (KMC) ver. " << KMC::CfgConsts::kmc_ver << " (" << KMC::CfgConsts::kmc_date << ")\n"
		<< "Usage:\n kmc [options] <input_file_name> <output_file_name> <working_directory>\n"
		<< " kmc [options] <@input_file_names> <output_file_name> <working_directory>\n"
		<< "Parameters:\n"
		<< "  input_file_name - single file in specified (-f switch) format (gziped or not)\n"
		<< "  @input_file_names - file name with list of input files in specified (-f switch) format (gziped or not)\n"
		<< "Options:\n"
		<< "  -v - verbose mode (shows all parameter settings); default: false\n"
		<< "  -k<len> - k-mer length (k from " << KMC::CfgConsts::min_k << " to " << KMC::CfgConsts::max_k << "; default: 25)\n"
		<< "  -m<size> - max amount of RAM in GB (from 1 to 1024); default: 12\n"
		<< "  -sm - use strict memory mode (memory limit from -m<n> switch will not be exceeded)\n"
		<< "  -hc - count homopolymer compressed k-mers (approximate and experimental)\n"
		<< "  -p<par> - signature length (5, 6, 7, 8, 9, 10, 11); default: 9\n"
		<< "  -f<a/q/m/bam/kmc> - input in FASTA format (-fa), FASTQ format (-fq), multi FASTA (-fm) or BAM (-fbam) or KMC(-fkmc); default: FASTQ\n"
		<< "  -ci<value> - exclude k-mers occurring less than <value> times (default: 2)\n"
		<< "  -cs<value> - maximal value of a counter (default: 255)\n"
		<< "  -cx<value> - exclude k-mers occurring more of than <value> times (default: 1e9)\n"
		<< "  -b - turn off transformation of k-mers into canonical form\n"
		<< "  -r - turn on RAM-only mode \n"
		<< "  -n<value> - number of bins \n"
		<< "  -t<value> - total number of threads (default: no. of CPU cores)\n"
		<< "  -sf<value> - number of FASTQ reading threads\n"
		<< "  -sp<value> - number of splitting threads\n"
		<< "  -sr<value> - number of threads for 2nd stage\n"
		<< "  -j<file_name> - file name with execution summary in JSON format\n"
		<< "  -w - without output\n"
		<< "  -o<kmc/kff> - output in KMC of KFF format; default: KMC\n"
		<< "  -hp - hide percentage progress (default: false)\n"
		<< "  -e - only estimate histogram of k-mers occurrences instead of exact k-mer counting\n"
		<< "  --opt-out-size - optimize output database size (may increase running time)\n"
		<< "Example:\n"
		<< "kmc -k27 -m24 NA19238.fastq NA.res /data/kmc_tmp_dir/\n"
		<< "kmc -k27 -m24 @files.lst NA.res /data/kmc_tmp_dir/\n";
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

bool CanCreateFile(const std::string& path)
{
	FILE* f = fopen(path.c_str(), "wb");
	if (!f)
		return false;
	fclose(f);
	remove(path.c_str());
	return true;
}

bool CanCreateFileInPath(const std::string& path)
{
	static const std::string name = "kmc_test.bin"; //Some random name
	if (path.back() == '\\' || path.back() == '/')
		return CanCreateFile(path + name);
	else
		return CanCreateFile(path + '/' + name);
}

void fill_temporary_kmc_databases_names(Params& params)
{
	for (uint32_t tmp_database_id = 0; tmp_database_id < params.mkmcParams.inputFiles.size(); ++tmp_database_id)
	{
		std::ostringstream sstreamKMC, sstreamTools;
		sstreamKMC << params.stage1ParamsTemplate.GetTmpPath();
		sstreamTools << params.stage1ParamsTemplate.GetTmpPath();
		if (params.stage1ParamsTemplate.GetTmpPath().back() != '/' && params.stage1ParamsTemplate.GetTmpPath().back() != '\\')
		{
			sstreamKMC << "/";
			sstreamTools << "/";
		}
		sstreamKMC << "kmc_db_" << std::setfill('0') << std::setw(5) << tmp_database_id;
		sstreamTools << "tools_db_" << std::setfill('0') << std::setw(5) << tmp_database_id;

		params.mkmcParams.kmcOutputFiles.push_back(sstreamKMC.str());
		params.mkmcParams.toolsOutputFiles.push_back(sstreamTools.str());
	}
}

//----------------------------------------------------------------------------------
// Parse the parameters
bool parse_parameters(int argc, char* argv[], Params& params)
{
	KMC::Stage1Params& stage1Params = params.stage1ParamsTemplate;
	KMC::Stage2Params& stage2Params = params.stage2ParamsTemplate;
	MKMCParams& mkmcParams = params.mkmcParams;
	int i;

	bool was_sm = false;
	bool was_r = false;

	bool was_e = false;
	bool was_opt_out_size = false;
	if (argc < 4)
		return false;

	for (i = 1; i < argc; ++i)
	{
		if (argv[i][0] != '-')
			break;
		// Number of threads
		if (strncmp(argv[i], "-t", 2) == 0)
		{
			mkmcParams.nThreads = atoi(&argv[i][2]);
		}
		// Number of parallel KMC runs
		else if (strncmp(argv[i], "-wrk", 4) == 0)
		{
			mkmcParams.nKMCWorkers = atoi(&argv[i][4]);
		}
		// k-mer length
		else if (strncmp(argv[i], "-k", 2) == 0)
			stage1Params.SetKmerLen(atoi(&argv[i][2]));
		// Memory limit
		else if (strncmp(argv[i], "-m", 2) == 0)
		{
			mkmcParams.maxRamGB = atoi(&argv[i][2]);
		}
		// Minimum counter threshold
		else if (strncmp(argv[i], "-ci", 3) == 0)
			stage2Params.SetCutoffMin(atoi(&argv[i][3]));
		// Maximum counter threshold
		else if (strncmp(argv[i], "-cx", 3) == 0)
			stage2Params.SetCutoffMax(atoll(&argv[i][3]));
		// Maximal counter value
		else if (strncmp(argv[i], "-cs", 3) == 0)
			stage2Params.SetCounterMax(atoll(&argv[i][3]));
		// Set p1
		else if (strncmp(argv[i], "-p", 2) == 0)
			stage1Params.SetSignatureLen(atoi(&argv[i][2]));
		//output type
		else if (strncmp(argv[i], "-o", 2) == 0)
		{
			if (strncmp(argv[i] + 2, "kff", 3) == 0)
				stage2Params.SetOutputFileType(KMC::OutputFileType::KFF);
			else if (strncmp(argv[i] + 2, "kmc", 3) == 0)
				stage2Params.SetOutputFileType(KMC::OutputFileType::KMC);
			else
			{
				std::cerr << "Error: unsupported output type: " << argv[i] << " (use -okff or -okmc)\n";
				exit(1);
			}
		}
		// FASTA input files
		else if (strncmp(argv[i], "-fa", 3) == 0)
			stage1Params.SetInputFileType(KMC::InputFileType::FASTA);
		// FASTQ input files
		else if (strncmp(argv[i], "-fq", 3) == 0)
			stage1Params.SetInputFileType(KMC::InputFileType::FASTQ);
		else if (strncmp(argv[i], "-fm", 3) == 0)
			stage1Params.SetInputFileType(KMC::InputFileType::MULTILINE_FASTA);
		else if (strncmp(argv[i], "-fbam", 5) == 0)
			stage1Params.SetInputFileType(KMC::InputFileType::BAM);
		else if (strncmp(argv[i], "-fkmc", 5) == 0)
			stage1Params.SetInputFileType(KMC::InputFileType::KMC);
#ifdef DEVELOP_MODE
		else if (strncmp(argv[i], "-vl", 3) == 0)
			stage1Params.SetDevelopVerbose(true);
#endif
		else if (strncmp(argv[i], "-v", 2) == 0)
		{
			static KMC::CerrVerboseLogger logger;
			stage1Params.SetVerboseLogger(&logger);
		}
		else if (strncmp(argv[i], "-sm", 3) == 0 && strlen(argv[i]) == 3)
		{
			was_sm = true;
			stage2Params.SetStrictMemoryMode(true);
		}
		else if (strncmp(argv[i], "-hc", 3) == 0 && strlen(argv[i]) == 3)
			stage1Params.SetHomopolymerCompressed(true);
		else if (strncmp(argv[i], "-r", 2) == 0)
		{
			stage1Params.SetRamOnlyMode(true);
			was_r = true;
		}
		else if (strncmp(argv[i], "-b", 2) == 0)
			stage1Params.SetCanonicalKmers(false);
		// Number of reading threads
		else if (strncmp(argv[i], "-sf", 3) == 0)
			stage1Params.SetNReaders(atoi(&argv[i][3]));
		// Number of splitting threads
		else if (strncmp(argv[i], "-sp", 3) == 0)
			stage1Params.SetNSplitters(atoi(&argv[i][3]));
		// Number of internal threads per 2nd stage
		else if (strncmp(argv[i], "-sr", 3) == 0)
			stage2Params.SetNThreads(atoi(&argv[i][3]));
		else if (strncmp(argv[i], "-n", 2) == 0)
			stage1Params.SetNBins(atoi(&argv[i][2]));
		/*else if (strncmp(argv[i], "-j", 2) == 0)
		{
			cliParams.jsonSummaryFileName = &argv[i][2];
			if (cliParams.jsonSummaryFileName == "")
				cerr << "Warning: file name for json summary file missed (-j switch)\n";
		}
		*/
		else if (strncmp(argv[i], "-e", 2) == 0)
		{
			was_e = true;
			stage1Params.SetEstimateHistogramCfg(KMC::EstimateHistogramCfg::ONLY_ESTIMATE);
		}
		else if (strcmp(argv[i], "--opt-out-size") == 0)
		{
			was_opt_out_size = true;
			if (stage1Params.GetEstimateHistogramCfg() != KMC::EstimateHistogramCfg::ONLY_ESTIMATE) //ONLY_ESTIMATE has priority over estimate and count
				stage1Params.SetEstimateHistogramCfg(KMC::EstimateHistogramCfg::ESTIMATE_AND_COUNT_KMERS);
		}
		else if (strncmp(argv[i], "-w", 2) == 0)
			stage2Params.SetWithoutOutput(true);

		if (strncmp(argv[i], "-smso", 5) == 0)
			stage2Params.SetStrictMemoryNSortingThreadsPerSorters(atoi(&argv[i][5]));

		if (strncmp(argv[i], "-smun", 5) == 0)
			stage2Params.SetStrictMemoryNUncompactors(atoi(&argv[i][5]));
		if (strncmp(argv[i], "-smme", 5) == 0)
			stage2Params.SetStrictMemoryNMergers(atoi(&argv[i][5]));

		if (strncmp(argv[i], "-dmp", 4) == 0)
			mkmcParams.dumpStepSize = atoll(&argv[i][4]);
		else if (strncmp(argv[i], "-hdd", 4) == 0 && strlen(argv[i]) == 4)
			mkmcParams.dumpToFile = true;
		else if (strncmp(argv[i], "-thr", 4) == 0)
			mkmcParams.minKmersPresenceThreshold = atof(&argv[i][4]);
	}

	if (argc - i < 3)
		return false;

	std::string input_file_name = std::string(argv[i++]);

	mkmcParams.outputFile = argv[i++];

	stage1Params.SetTmpPath(argv[i++]);

	std::vector<std::string> input_file_names;
	if (input_file_name[0] != '@')
		input_file_names.push_back(input_file_name);
	else
	{
		std::ifstream in(input_file_name.c_str() + 1);
		if (!in.good())
		{
			std::cerr << "Error: No " << input_file_name.c_str() + 1 << " file\n";
			return false;
		}

		std::string s;
		while (std::getline(in, s))
			if (s != "")
				input_file_names.push_back(s);

		in.close();
	}
	mkmcParams.inputFiles.swap(input_file_names);

	fill_temporary_kmc_databases_names(params);

	//Validate and resolve conflicts in parameters
	if (was_e && was_opt_out_size)
	{
		std::cerr << "Warning: --opt-out-size is ignored because -e was used\n";
	}

	if (was_sm && was_r)
	{
		std::cerr << "Error: -sm can not be used with -r\n";
		return false;
	}

	//Check if output files may be created and if it is possible to create file in specified tmp location
	if (!stage2Params.GetWithoutOutput())
	{
		std::string pre_file_name = stage2Params.GetOutputFileName() + ".kmc_pre";
		std::string suff_file_name = stage2Params.GetOutputFileName() + ".kmc_suf";
		if (!CanCreateFile(pre_file_name))
		{
			std::cerr << "Error: Cannot create file: " << pre_file_name << "\n";
			return false;
		}
		if (!CanCreateFile(suff_file_name))
		{
			std::cerr << "Error: Cannot create file: " << suff_file_name << "\n";
			return false;
		}
	}
	if (!CanCreateFileInPath(stage1Params.GetTmpPath()))
	{
		std::cerr << "Error: Cannot create file in specified working directory: " << stage1Params.GetTmpPath() << "\n";
		return false;
	}
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
		params.setKMCParams();

		std::cout << "Starting k-mer counting...\n";
		KMCRunner kmcRunner(params);
		kmcRunner.runKMCParallel();

		std::cout << "Starting converting k-mer databases...\n";
		KMCToolsRunner kmcToolsRunner(params);
		kmcToolsRunner.runKMCToolsParallel();

		Dump dump(params);
		if (params.mkmcParams.dumpToFile)
		{
			std::cout << "Starting dumping to file " << params.mkmcParams.outputFile << "...\n";
			dump.dumpToFile();
		}
		else
		{
			std::cout << "Starting dumping to stdout...\n";
			dump.dumpToStd();
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << e.what() << '\n';
	}
}