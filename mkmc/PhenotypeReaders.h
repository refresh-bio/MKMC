#pragma once

#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <cstdint>


struct Params;

namespace PhenotypeReaderHelpers
{
	size_t get_no_samples(const Params& params);
}


template<typename Phenotype_T>
class PhenotypeReader
{
	std::string fileName;
	std::vector<Phenotype_T> phenotype;

protected:
	const Params& params;

public:
	PhenotypeReader(const Params& params) :
		params(params)
	{}

	void setFileName(const std::string& _fileName) { fileName = _fileName; }
	const std::string& getFileName() const { return fileName; }

	bool readPhenotype();
	const std::vector<Phenotype_T>& getPhenotype() const { return phenotype; }
};



class DifferentialAnalysisPhenotypeReader : public PhenotypeReader<std::string>
{
	std::map<std::string, uint32_t> mapToInt;
	std::vector<std::string> mapToClass;

	std::vector<uint32_t> mappedPhenotype;

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
	DifferentialAnalysisPhenotypeReader(const Params& params) :
		PhenotypeReader<std::string>(params)
	{}

	bool mapPhenotypeToInts();

	const std::vector<uint32_t>& getMappedPhenotype() const { return mappedPhenotype; }
	size_t getClassesNumber() const { return mapToInt.size(); };
};



template<typename Phenotype_T>
bool PhenotypeReader<Phenotype_T>::readPhenotype()
{
	std::ifstream phenotypeFile(fileName);
	if (!phenotypeFile.is_open())
	{
		std::cerr << "Error: cannot open " << fileName << "." << std::endl;
		return false;
	}

	Phenotype_T value;
	while (phenotypeFile >> value)
		phenotype.push_back(value);

	if (phenotypeFile.fail() && !phenotypeFile.eof())
	{
		std::cerr << "Error: wrong value in a file " << fileName << "." << std::endl;
		return false;
	}

	auto no_samples = PhenotypeReaderHelpers::get_no_samples(params);
	if (phenotype.size() != no_samples)
	{
		std::cerr << "Error: number of a phenotype values in a file " << fileName << " (" << phenotype.size() << ") is different than number of samples (" << no_samples << ")." << std::endl;
		return false;
	}

	return true;
}
