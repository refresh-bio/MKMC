#pragma once

#include <mutex>
#include <vector>
#include <cstdint>
#include <algorithm>
#include "parameters.h"



class TasksPool
{
	uint32_t nTasks;
	uint32_t nextTask;
	std::mutex taskAvailableMutex;

public:
	TasksPool(const uint32_t nTasks) :
		nTasks(nTasks),
		nextTask(0)
	{}

	bool getTask(uint32_t& task);
};

class KMCRunner
{
	const Params& params;

	TasksPool tasksPool;
	void operator()();

public:
	KMCRunner(const Params& params) :
		params(params),
		tasksPool(static_cast<uint32_t>(params.mkmcParams.inputFiles.size()))
	{}

	void runKMCParallel();
};
