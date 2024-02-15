#include <thread>
#include <vector>
#include <cstdint>
#include <filesystem>
#include <limits>
#include "KMCRunner.h"
#include "Logger.h"



void KMCRunner::operator()()
{
	uint32_t task = std::numeric_limits<uint32_t>::max();
	while (tasksPool.getTask(task))
	{
		std::string inputFile = params.mkmcParams.inputFiles[task];
		std::string outputFile = params.mkmcParams.kmcOutputFiles[task];
		std::string tmpDir = params.mkmcParams.kmcTmpDirs[task];

		Logger::Inst().Log("Start k-mer counting for " + inputFile);
		KMC::Runner runner;

		// Fill missing, per counting, KMC parameters
		KMC::Stage1Params stage1Params = params.stage1Params;
		stage1Params.SetInputFiles({ inputFile });
		stage1Params.SetTmpPath(tmpDir);
		runner.RunStage1(stage1Params);

		KMC::Stage2Params stage2Params = params.stage2Params;
		stage2Params.SetOutputFileName(outputFile);
		runner.RunStage2(stage2Params);
		Logger::Inst().Log("k-mer counting for " + inputFile + " done.");
	}
}

void KMCRunner::runKMCParallel()
{
	if (!params.stage1Params.GetRamOnlyMode())
		for (const auto& dirPath : params.mkmcParams.kmcTmpDirs)
			std::filesystem::create_directory(dirPath);

	std::vector<std::thread> threads(params.mkmcParams.nKMCWorkers);
	for (uint32_t i_thred = 0; i_thred < params.mkmcParams.nKMCWorkers; ++i_thred)
	{
		threads[i_thred] = std::thread([this] { (*this)(); });
	}

	for (std::thread& thread : threads)
	{
		thread.join();
	}

	if (!params.stage1Params.GetRamOnlyMode())
		for (const auto& dirPath : params.mkmcParams.kmcTmpDirs)
			std::filesystem::remove(dirPath);
}
