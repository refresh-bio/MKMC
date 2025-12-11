#include "SamplesFileReader.h"
#include "Logger.h"
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
		Logger::Inst().Log("Error: a sample " + sampleName + " is given multiple times in " + mkmcParams.inputFileName + ".");
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
		Logger::Inst().Log("Error: cannot open " + fileName + " (" + mkmcParams.inputFileName + ", line " + std::to_string(lineNo) + ").");
		return false;
	}
	return true;
}



bool SamplesFileReader::readSamples(std::vector<Sample>& oSamples, bool& warningPrinted)
{
	std::ifstream in(mkmcParams.inputFileName);
	if (!in.good())
	{
		Logger::Inst().Log("Error: no " + mkmcParams.inputFileName + " file.");
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
		Logger::Inst().Log("Error: no samples specified in a " + mkmcParams.inputFileName + " file.");
		return false;
	}


	if (singleWordLines)
	{
		Logger::Inst().Log("Warning: some of lines in input file " + mkmcParams.inputFileName + " contain just one word, they will be treated both as samples names and files names.", 1);
		warningPrinted = true;
	}

	oSamples = samples;
	return true;
}
