#pragma once

#include <mutex>
#include <vector>
#include <cstdint>
#include <algorithm>
#include <string>
#define NOMINMAX
#include "progress_bar.hpp"
#include "parameters.h"
#include "TasksPool.h"


class KMCRunner
{
	struct TaskData
	{
		std::vector<std::string> inputFiles;
		std::string outputFile;
		std::string tmpDir;
		KMC::InputFileType inputFileType = KMC::InputFileType::FASTQ;
	};
	std::vector<TaskData> tasksData;

	const Params& params;

	TasksPool<TaskData> tasksPool;

	ProgressBar progress_bar;

	void operator()();

public:
	KMCRunner(const Params& params) :
		params(params),
		tasksPool(tasksData),
		progress_bar((params.mkmcParams.samples.size() + params.filterParams.filterKmersSequences) * 200, "k-mer counting", std::cerr, params.mkmcParams.verbosity_level == 0)
	{
		const MKMCParams& mkmcParams = params.mkmcParams;
		tasksData.reserve(mkmcParams.samples.size());
		for (uint32_t i = 0; i < mkmcParams.samples.size(); ++i)
		{
			tasksData.push_back(TaskData{ mkmcParams.samples[i].inputFiles, mkmcParams.kmcOutputFiles[i], mkmcParams.kmcTmpDirs[i], params.mkmcParams.inputFileType });
		}

		if (params.filterParams.filterKmersSequences)
		{
			tasksData.push_back(TaskData{ { params.mutableParams.kmersSequencesToFilterOut }, params.filterParams.kmersSequencesToFilterOutDB, mkmcParams.tmpPath, KMC::InputFileType::MULTILINE_FASTA });
		}
	}

	void runKMCParallel();
};
