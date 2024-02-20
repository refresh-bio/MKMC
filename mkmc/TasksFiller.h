#pragma once

#include "parameters.h"
#include <string>
#include <vector>



class TasksFiller
{
	const MKMCParams& mkmcParams;
	const std::string inputFileName;

	std::vector<std::string> samplesNames;
	std::vector<std::vector<std::string>> inputFilesPerSample;

	bool parseLine(const std::string& line, uint32_t lineNo);
	bool canOpenFile(const std::string& fileName, uint32_t lineNo);

public:
	TasksFiller(const MKMCParams& mkmcParams, const std::string& inputFileName) :
		mkmcParams(mkmcParams),
		inputFileName(inputFileName)
	{}

	bool readSamples(std::vector<std::string>& oSamplesNames, std::vector<std::vector<std::string>>& oInputFilesPerSample);
};