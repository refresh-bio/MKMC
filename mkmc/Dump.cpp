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
	uint32_t signatureLen = samples.front().GetSignatureLen();
	auto signatureSelectionScheme = samples.front().GetSignatureSelectionScheme();

	for (const auto& sample : samples)
	{
		if (k != sample.GetK())
			return false;
		if (signatureLen != sample.GetSignatureLen())
			return false;
		if (signatureSelectionScheme != sample.GetSignatureSelectionScheme())
			return false;
	}
	return true;
}



void Dump::openDatabases(std::vector<KMCFileWrapper>& samples, uint32_t binId)
{
	for (const std::string& fileName : params.mkmcParams.kmcOutputFiles)
	{
		samples.emplace_back(fileName, binId);
	}

	if (!allKAreSame(samples))
	{
		std::cerr << "Error: KMC databases are not consistent." << std::endl;
		exit(1);
	}
}



void Dump::dumpToFileParallel()
{
	tasksData.reserve(params.stage1Params.GetNBins());
	for (uint32_t i = 0; i < params.stage1Params.GetNBins(); ++i)
	{
		tasksData.push_back(TaskData{ i });
	}

	std::vector<std::thread> threads(params.mkmcParams.nDumpThreads);
	for (uint32_t i_thred = 0; i_thred < params.mkmcParams.nDumpThreads; ++i_thred)
	{
		threads[i_thred] = std::thread([this] { (*this)(); });
	}

	for (std::thread& thread : threads)
	{
		thread.join();
	}
}



void Dump::operator()()
{
	TaskData taskData;
	if (params.mkmcParams.outputFileType == OutputFileType::Matrix)
		while (tasksPool.getTask(taskData))
		{
			dumpToFile<MatrixFileGenerator>(params.mkmcParams.outputFiles[taskData.binId], taskData.binId);
		}
	else if (params.mkmcParams.outputFileType == OutputFileType::FASTA)
		while (tasksPool.getTask(taskData))
		{
			dumpToFile<FASTAFileGenerator>(params.mkmcParams.outputFiles[taskData.binId], taskData.binId);
		}
	else
		assert(false);
}