#pragma once

#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <sstream>



struct Params;
class DifferentialAnalysisPhenotypeParser
{
	const Params& params;

	std::map<std::string, uint32_t> mapToInt;
	std::vector<std::string> mapToClass;

	std::vector<std::string> samplesClasses;
	std::vector<uint32_t> mappedClasses;

	bool isNatural(const std::string& str)
	{
		return std::find_if(str.begin(), str.end(), [](unsigned char c) { return !std::isdigit(c); }) == str.end();
	}

	uint32_t strToNatural(const std::string& str)
	{
		uint32_t res;
		std::istringstream sstream(str);
		sstream >> res;
		return res;
	}

	bool verifyClassesSense();

public:
	DifferentialAnalysisPhenotypeParser(const Params& params) :
		params(params)
	{}

	void readClasses();
	void mapClassesToInt();

	const std::vector<uint32_t>& getMappedClasses() const { return mappedClasses; }
};
