#include "TasksFiller.h"
#include <fstream>
#include <sstream>
#include <string>



bool TasksFiller::parseLine(const std::string& line, uint32_t lineNo)
{
	std::istringstream lineStream(line);
	std::string sampleName, fileName;
	lineStream >> sampleName;
	if (!lineStream) // nothing in line
		return true;

	samplesNames.push_back(sampleName);
	inputFilesPerSample.push_back(std::vector<std::string>());

	lineStream >> fileName;
	if (!lineStream) // sample name = file name
	{
		if (mkmcParams.verbosity_level > 0)
			std::cerr << "Warning: input file line " << lineNo << " contains just one word, it will be treated both as a sample and a file name." << std::endl;

		if (!canOpenFile(sampleName, lineNo))
			return false;

		inputFilesPerSample.back().push_back(sampleName);
	}
	else
	{
		do
		{
			if (!canOpenFile(fileName, lineNo))
				return false;

			inputFilesPerSample.back().push_back(fileName);

			lineStream >> fileName;
		} while (lineStream);
	}

	return true;
}



bool TasksFiller::canOpenFile(const std::string& fileName, uint32_t lineNo)
{
	std::ifstream inFile(fileName);
	if (!inFile.is_open())
	{
		std::cerr << "Error: Cannot open " << fileName << " (" << mkmcParams.inputFileName << ", line " << lineNo << ").\n" << std::endl;
		return false;
	}
	return true;
}



bool TasksFiller::readSamples(std::vector<std::string>& oSamplesNames, std::vector<std::vector<std::string>>& oInputFilesPerSample)
{
	std::ifstream in(mkmcParams.inputFileName);
	if (!in.good())
	{
		std::cerr << "Error: No " << mkmcParams.inputFileName << " file.\n" << std::endl;
		return false;
	}

	std::string line;
	uint32_t lineNo = 1;
	while (std::getline(in, line))
	{
		if (!parseLine(line, lineNo))
			return false;
		++lineNo;
	}

	oSamplesNames.swap(samplesNames);
	oInputFilesPerSample.swap(inputFilesPerSample);
	return true;
}
