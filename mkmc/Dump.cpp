#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstdint>
#include <limits>
#include <memory>
#include <sstream>
#include <algorithm>
#include <numeric>
#include "../kmc/kmc_api/kmc_file.h"
#include "../kmc/kmc_dump/nc_utils.h"
#include "Dump.h"
#include "Filter.h"
#include "FileGenerators.h"



bool Dump::allKAreSame(const std::vector<KMCFileWrapper>& samples)
{
	if (samples.empty())
		return true;
	uint32_t k = samples.front().GetK();
	for (const auto& sample : samples)
		if (k != sample.GetK())
			return false;
	return true;
}

void Dump::openDatabases(std::vector<KMCFileWrapper>& samples)
{
	
	for (const std::string& fileName : params.mkmcParams.toolsOutputFiles)
	{
		samples.emplace_back(fileName);
	}

	if (!allKAreSame(samples))
	{
		std::cerr << "Error: each database should have the same k." << std::endl;
		exit(1);
	}
}
