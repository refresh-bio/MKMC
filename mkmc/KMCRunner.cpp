#include <thread>
#include <vector>
#include "KMCRunner.h"



void KMCRunner::operator()()
{
	std::string inputFile, outputFile;
	while (tasksPool.getTask(inputFile, outputFile))
	{
		KMC::Runner runner;

		KMC::Stage1Params stage1Params = params.stage1ParamsTemplate;
		stage1Params.SetInputFiles({ inputFile });
		auto stage1Result = runner.RunStage1(stage1Params);

		KMC::Stage2Params stage2Params = params.stage2ParamsTemplate;
		stage2Params.SetOutputFileName(outputFile);
		auto stage2Result = runner.RunStage2(stage2Params);
	}
}

void KMCRunner::runKMCParallel()
{
	std::vector<std::thread> threads(params.mkmcParams.nKMCWorkers);
	for (uint32_t i_thred = 0; i_thred < params.mkmcParams.nKMCWorkers; ++i_thred)
	{
		threads[i_thred] = std::thread([this] { (*this)(); });
	}

	for (std::thread& thread : threads)
	{
		thread.join();
	}
}

bool TasksPool::getTask(std::string& inputFile, std::string& outputFile) {
	std::unique_lock<std::mutex> lck(taskAvailableMutex);
	assert(nextTask <= inputFiles.size());
	if (nextTask == inputFiles.size())
	{
		return false;
	}
	inputFile = inputFiles[nextTask];
	outputFile = outputFiles[nextTask];
	++nextTask;
	return true;
}
