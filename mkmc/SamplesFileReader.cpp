#include "SamplesFileReader.h"
#include <fstream>
#include <sstream>
#include <string>



bool SamplesFileReader::parseLine(const std::string& line, uint32_t lineNo, bool& singleWordLine)
{
	std::istringstream lineStream(line);
	std::string sampleName, fileName;
	lineStream >> sampleName;

	if (!lineStream) // nothing in line
		return true;

	if (samplesNames.find(sampleName) != samplesNames.end())
	{
		std::cerr << "Error: Sample " << sampleName << " is given multiple times in " << mkmcParams.inputFileName << "." << std::endl;
		return false;
	}
	samplesNames.insert(sampleName);

	samples.push_back(Sample{ sampleName, std::vector<std::string>() });

	lineStream >> fileName;
	if (!lineStream) // sample name = file name
	{
		singleWordLine = true;

		if (!canOpenFile(sampleName, lineNo))
			return false;

		samples.back().inputFiles.push_back(sampleName);
	}
	else
	{
		do
		{
			if (!canOpenFile(fileName, lineNo))
				return false;

			samples.back().inputFiles.push_back(fileName);

			lineStream >> fileName;
		} while (lineStream);
	}

	return true;
}



bool SamplesFileReader::canOpenFile(const std::string& fileName, uint32_t lineNo)
{
	std::ifstream inFile(fileName);
	if (!inFile.is_open())
	{
		std::cerr << "Error: Cannot open " << fileName << " (" << mkmcParams.inputFileName << ", line " << lineNo << ")." << std::endl;
		return false;
	}
	return true;
}



bool SamplesFileReader::readSamples(std::vector<Sample>& oSamples, bool& warningPrinted)
{
	std::ifstream in(mkmcParams.inputFileName);
	if (!in.good())
	{
		std::cerr << "Error: No " << mkmcParams.inputFileName << " file." << std::endl;
		return false;
	}

	std::string line;
	uint32_t lineNo = 1;
	bool singleWordLines = false;
	while (std::getline(in, line))
	{
		if (!parseLine(line, lineNo, singleWordLines))
			return false;
		++lineNo;
	}

	if (samples.empty())
	{
		std::cerr << "Error: No samples specified in a " << mkmcParams.inputFileName << " file." << std::endl;
		return false;
	}


	if (singleWordLines && mkmcParams.verbosity_level > 0)
	{
		std::cerr << "Warning: some of lines in input file " << mkmcParams.inputFileName << " contain just one word, they will be treated both as samples names and files names." << std::endl;
		warningPrinted = true;
	}

	oSamples = samples;
	return true;
}
