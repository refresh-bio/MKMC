#pragma once

#include <mutex>
#include <vector>
#include <cstdint>
#include <algorithm>
#include <string>
#include "parameters.h"
#include "TasksPool.h"



class KMCRunner
{
	struct TaskData
	{
		std::vector<std::string> inputFiles;
		std::string outputFile;
		std::string tmpDir;
	};
	std::vector<TaskData> tasksData;

	const Params& params;

	TasksPool<TaskData> tasksPool;
	void operator()();

public:
	KMCRunner(const Params& params) :
		params(params),
		tasksPool(tasksData)
	{
		const MKMCParams& mkmcParams = params.mkmcParams;
		tasksData.reserve(mkmcParams.inputFilesPerSample.size());
		for (uint32_t i = 0; i < mkmcParams.inputFilesPerSample.size(); ++i)
		{
			tasksData.push_back(TaskData{ mkmcParams.inputFilesPerSample[i], mkmcParams.kmcOutputFiles[i], mkmcParams.kmcTmpDirs[i]});
		}
	}

	void runKMCParallel();
};
