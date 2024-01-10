#pragma once

#include <mutex>
#include "parameters.h"



class TasksPool
{
	const Params& params;
	uint32_t nextTask;
	std::mutex taskAvailableMutex;

public:
	TasksPool(const Params& params) :
		params(params),
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
		tasksPool(params)
	{}

	void runKMCParallel();
};
