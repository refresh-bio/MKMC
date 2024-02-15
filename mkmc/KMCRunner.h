#pragma once

#include <mutex>
#include <vector>
#include <cstdint>
#include <algorithm>
#include "parameters.h"
#include "TasksPool.h"



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
