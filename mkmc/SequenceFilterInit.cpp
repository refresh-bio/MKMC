#include "SequenceFilterInit.h"



bool SequenceFilterInit::getNotEmptyLine(std::istream& stream, std::string& outLine)
{
	std::string line;
	while (!stream.eof())
	{
		std::getline(stream, line);
		if (!line.empty())
		{
			outLine = line;
			return true;
		}
	}
	return false;
}



bool SequenceFilterInit::isFasta()
{
	std::ifstream stream(params.filterParams.inputKmersSequencesToFilterOut);
	if (!stream.is_open())
	{
		std::cerr << "Error: cannot open " << params.filterParams.inputKmersSequencesToFilterOut << "." << std::endl;
		exit(1);
	}
	std::string line;
	std::string kmer;

	if (!getNotEmptyLine(stream, line))
	{
		std::cerr << "Error: format of a file " << params.filterParams.inputKmersSequencesToFilterOut << "has not proper format; it must be one of: FASTA or a sequence of k-mers." << std::endl;
		exit(1);
	}

	if (line.front() == '>')
	{
		if (!getNotEmptyLine(stream, kmer))
		{
			std::cerr << "Error: k-mers to be filtered out have inproper length, it must equal " << params.stage1Params.GetKmerLen() << "." << std::endl;
			exit(1);
		}

		if (kmer.length() != params.stage1Params.GetKmerLen())
		{
			std::cerr << "Error: k-mers to be filtered out have inproper length, it must equal " << params.stage1Params.GetKmerLen() << "." << std::endl;
			exit(1);
		}
		return true;
	}
	else
	{
		std::istringstream sstream(line);
		sstream >> kmer;
		if (kmer.length() != params.stage1Params.GetKmerLen())
		{
			std::cerr << "Error: k-mers to be filtered out have inproper length, it must equal " << params.stage1Params.GetKmerLen() << "." << std::endl;
			exit(1);
		}
		return false;
	}
}



void SequenceFilterInit::convertTxtToFasta()
{
	std::ofstream filterFasta(params.mutableParams.kmersSequencesToFilterOut);
	std::ifstream filterTxt(params.filterParams.inputKmersSequencesToFilterOut);

	if (!filterFasta.is_open())
	{
		std::cerr << "Error: cannot create temporary file to perform filtering k-mers sequences out." << std::endl;
		exit(1);
	}

	if (!filterTxt.is_open())
	{
		std::cerr << "Error: cannot open " << params.filterParams.inputKmersSequencesToFilterOut << "." << std::endl;
		exit(1);
	}

	std::string kmer;
	while (!filterTxt.eof())
	{
		filterTxt >> kmer;
		if (kmer.empty())
			break;

		if (kmer.length() != params.stage1Params.GetKmerLen())
		{
			std::cerr << "Error: k-mers to be filtered out have inproper length, it must equal " << params.stage1Params.GetKmerLen() << "." << std::endl;
			exit(1);
		}
		filterFasta << ">\n" << kmer << '\n';
		kmer.clear();
	}
}



void SequenceFilterInit::prepareKmersSequencesToFilter()
{
	if (!isFasta())
	{
		std::cerr << "Starting preparing k-mers to filter out..." << std::endl;

		sequence_filter_init.startTimer();
		convertTxtToFasta();
		sequence_filter_init.stopTimer();

		params.mutableParams.createdFastaFile = true;
	}
	else
	{
		params.mutableParams.kmersSequencesToFilterOut = params.filterParams.inputKmersSequencesToFilterOut;
	}
}
