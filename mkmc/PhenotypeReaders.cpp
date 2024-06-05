#include "PhenotypeReaders.h"
#include "parameters.h"
#include <string>



bool DifferentialAnalysisPhenotypeReader::verifyClassesSense()
{
	bool contains0 = false;
	uint32_t maximal = 0;
	for (const auto& c : mapToClass)
	{
		if (!isNatural(c))
			return true; // does not contains just numbers

		uint32_t n = strToNatural(c);
		if (n == 0)
			contains0 = true;
		if (n > maximal)
			maximal = n;
	}
	// if here, classes are natural numbers; perform remapping for convenience
	if (contains0 && maximal == mapToClass.size() - 1) // mapToClass contains unique values: two of them are 0 and n-1, hence rest of them are from a range of <1, n-2>
	{
		// For simplicity: if classes are natural numbers, map them basing on classes values
		mapToInt.clear();
		mapToClass.clear();

		for (uint32_t i = 0; i < maximal; ++i)
		{
			std::string classStr = std::to_string(i);
			mapToInt[classStr] = i;
			mapToClass.push_back(classStr);
		}
		for (uint32_t i = 0; i < maximal; ++i)
		{
			mappedPhenotype[i] = mapToInt[getPhenotype()[i]];
		}
		return true;
	}

	return false;
}



void DifferentialAnalysisPhenotypeReader::mapPhenotypeToInts()
{
	mappedPhenotype.reserve(getPhenotype().size());

	for (const auto& c : getPhenotype())
	{
		if (!mapToInt.contains(c))
		{
			mapToClass.push_back(c);
			mapToInt[c] = mapToClass.size() - 1;

			mappedPhenotype.push_back(mapToClass.size() - 1);
		}
		else
		{
			mappedPhenotype.push_back(mapToInt[c]);
		}
	}

	if (mapToInt.size() == 1)
	{
		std::cerr << "Error: a number of distinct classes in a file " << getFileName() << " must be greater than 1." << std::endl;
		exit(1);
	}

	if (!verifyClassesSense())
	{
		std::cerr << "Error: a file " << getFileName() << " must contain classes for all samples which must be natural numbers starting from 0 or text labels." << std::endl;
		exit(1);
	}
}
