#include <thread>
#include <vector>
#include <cstdint>
#include <filesystem>
#include <limits>
#include "KMCRunner.h"

class KMCPercentProgressObserver : public KMC::IPercentProgressObserver
{
	std::atomic<uint64_t> prev_val{};
	ProgressBar& progress_bar;
	void SetLabel(const std::string& label) override
	{
		//ignore labels
	}
	void ProgressChanged(int newValue) override
	{
		progress_bar += newValue - prev_val;
		prev_val = newValue;
	}
public:
	KMCPercentProgressObserver(ProgressBar& progress_bar) :
		progress_bar(progress_bar)
	{

	}
	void reset()
	{
		prev_val = 0;
	}
};

void KMCRunner::operator()()
{
	TaskData taskData;
	KMCPercentProgressObserver progress_observer(progress_bar);
	while (tasksPool.getTask(taskData))
	{
		std::string inputFiles;
		for (const auto& inputFile : taskData.inputFiles)
			inputFiles += " " + inputFile;

		KMC::Runner runner;

		// Fill missing, per counting, KMC parameters
		KMC::Stage1Params stage1Params = params.stage1Params;
		stage1Params.SetInputFiles(taskData.inputFiles);
		stage1Params.SetTmpPath(taskData.tmpDir);
		stage1Params.SetInputFileType(taskData.inputFileType);
		progress_observer.reset();
		stage1Params.SetPercentProgressObserver(&progress_observer);

		runner.RunStage1(stage1Params);

		KMC::Stage2Params stage2Params = params.stage2Params;
		stage2Params.SetOutputFileName(taskData.outputFile);
		progress_observer.reset();

		runner.RunStage2(stage2Params);
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
