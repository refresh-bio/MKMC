#pragma once

#include <mutex>
#include <vector>
#include <cstdint>
#include "parameters.h"



class TasksPool
{
	const std::vector<std::string>& inputFiles;
	const std::vector<std::string>& outputFiles;
	uint32_t nextTask;
	std::mutex taskAvailableMutex;

public:
	TasksPool(const std::vector<std::string>& inputFiles, const std::vector<std::string>& outputFiles) :
		inputFiles(inputFiles),
		outputFiles(outputFiles),
		nextTask(0)
	{}

	bool getTask(std::string& inputFile, std::string& outputFile);
};

class KMCRunner
{
	const Params& params;

	TasksPool tasksPool;
	void operator()();

public:
	KMCRunner(const Params& params) :
		params(params),
		tasksPool(params.mkmcParams.inputFiles, params.mkmcParams.kmcOutputFiles)
	{}

	void runKMCParallel();
};
