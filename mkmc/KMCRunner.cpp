#include <thread>
#include <vector>
#include <cstdint>
#include <filesystem>
#include <limits>
#include "KMCRunner.h"
#include "Logger.h"



void KMCRunner::operator()()
{
	TaskData taskData;
	while (tasksPool.getTask(taskData))
	{
		std::string inputFiles;
		for (const auto& inputFile : taskData.inputFiles)
			inputFiles += " " + inputFile;

		Logger::Inst().Log("Start k-mer counting for" + inputFiles);
		KMC::Runner runner;

		// Fill missing, per counting, KMC parameters
		KMC::Stage1Params stage1Params = params.stage1Params;
		stage1Params.SetInputFiles(taskData.inputFiles);
		stage1Params.SetTmpPath(taskData.tmpDir);

		runner.RunStage1(stage1Params);

		KMC::Stage2Params stage2Params = params.stage2Params;
		stage2Params.SetOutputFileName(taskData.outputFile);

		runner.RunStage2(stage2Params);
		Logger::Inst().Log("k-mer counting for" + inputFiles + " done.");
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
