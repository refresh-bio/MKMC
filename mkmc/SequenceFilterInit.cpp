#include "SequenceFilterInit.h"
#include "Logger.h"



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



bool SequenceFilterInit::isFastaOrMultiFasta()
{
	std::ifstream stream(params.filterParams.inputKmersSequencesToFilterOut);
	if (!stream.is_open())
	{
		Logger::Inst().Log("Error: cannot open " + params.filterParams.inputKmersSequencesToFilterOut + ".");
		exit(1);
	}
	std::string line;

	if (!getNotEmptyLine(stream, line))
	{
		Logger::Inst().Log("Error: format of a file " + params.filterParams.inputKmersSequencesToFilterOut + "is not proper; it must be one of: FASTA or a sequence of k-mers.");
		exit(1);
	}

	if (line.front() == '>') // FASTA or multi-FASTA
	{
		std::string firstKmerSeq;
		if (!getNotEmptyLine(stream, firstKmerSeq))
		{
			Logger::Inst().Log("Error: format of a file " + params.filterParams.inputKmersSequencesToFilterOut + "is not proper; it must be one of: FASTA or a sequence of k-mers.");
			exit(1);
		}
		return true;
	}
	else // raw k-mers file
	{
		std::istringstream sstream(line);
		std::string kmer;
		sstream >> kmer;
		if (kmer.length() != params.stage1Params.GetKmerLen())
		{
			Logger::Inst().Log("Error: k-mers to be filtered out have inproper length, it must equal " + std::to_string(params.stage1Params.GetKmerLen()) + ".");
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
		Logger::Inst().Log("Error: cannot create temporary file to perform filtering k-mers sequences out.");
		exit(1);
	}

	if (!filterTxt.is_open())
	{
		Logger::Inst().Log("Error: cannot open " + params.filterParams.inputKmersSequencesToFilterOut + ".");
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
			Logger::Inst().Log("Error: k-mers to be filtered out have inproper length, it must equal " + std::to_string(params.stage1Params.GetKmerLen()) + ".");
			exit(1);
		}
		filterFasta << ">\n" << kmer << '\n';
		kmer.clear();
	}
}



bool SequenceFilterInit::prepareKmersSequencesToFilter()
{
	if (!isFastaOrMultiFasta())
	{
		Logger::Inst().Log("Starting preparing k-mers for filtering...");

		sequence_filter_init.startTimer();
		convertTxtToFasta();
		sequence_filter_init.stopTimer();

		params.mutableParams.createdFastaFile = true;
		return true;
	}
	else
	{
		params.mutableParams.kmersSequencesToFilterOut = params.filterParams.inputKmersSequencesToFilterOut;
		return false;
	}
}
