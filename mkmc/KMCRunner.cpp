#include <thread>
#include <vector>
#include <cstdint>
#include <filesystem>
#include <limits>
#include "KMCRunner.h"



void KMCRunner::operator()()
{
	uint32_t task = std::numeric_limits<uint32_t>::max();
	while (tasksPool.getTask(task))
	{
		std::string inputFile = params.mkmcParams.inputFiles[task];
		std::string outputFile = params.mkmcParams.kmcOutputFiles[task];
		std::string tmpDir = params.mkmcParams.kmcTmpDirs[task];

		KMC::Runner runner;

		KMC::Stage1Params stage1Params = params.stage1ParamsTemplate;
		stage1Params.SetInputFiles({ inputFile });
		stage1Params.SetTmpPath(tmpDir);
		auto stage1Result = runner.RunStage1(stage1Params);

		KMC::Stage2Params stage2Params = params.stage2ParamsTemplate;
		stage2Params.SetOutputFileName(outputFile);
		auto stage2Result = runner.RunStage2(stage2Params);
	}
}

void KMCRunner::runKMCParallel()
{
	if (!params.stage1ParamsTemplate.GetRamOnlyMode())
		for (const auto& dirPath : params.mkmcParams.kmcTmpDirs)
			std::filesystem::create_directory(std::filesystem::path(dirPath));

	std::vector<std::thread> threads(params.mkmcParams.nKMCWorkers);
	for (uint32_t i_thred = 0; i_thred < params.mkmcParams.nKMCWorkers; ++i_thred)
	{
		threads[i_thred] = std::thread([this] { (*this)(); });
	}

	for (std::thread& thread : threads)
	{
		thread.join();
	}

	if (!params.stage1ParamsTemplate.GetRamOnlyMode())
		for (const auto& dirPath : params.mkmcParams.kmcTmpDirs)
			std::filesystem::remove(std::filesystem::path(dirPath));
}
