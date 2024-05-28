#pragma once

#include "parameters.h"
#include <string>
#include <vector>
#include <set>



class SamplesFileReader
{
	const MKMCParams& mkmcParams;

	std::vector<Sample> samples;
	std::set<std::string> samplesNames; // to easy search

	bool parseLine(const std::string& line, uint32_t lineNo, bool& singleWordLines);
	bool canOpenFile(const std::string& fileName, uint32_t lineNo);

public:
	SamplesFileReader(const MKMCParams& mkmcParams) :
		mkmcParams(mkmcParams)
	{}

	bool readSamples(std::vector<Sample>& oSamples);
};
