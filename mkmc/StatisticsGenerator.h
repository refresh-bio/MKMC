#pragma once

#include <fstream>
#include <vector>
#include <string>
#include <memory>
#include "parameters.h"
#include "TasksPool.h"
#include "refresh/statistics/lib/statistics.h"
#define NOMINMAX
#include "progress_bar.hpp"
#include "kmcdb/kmcdb.h"
#include "DumpWriter.h"
#include "StatisticsGatherers.h"


class StatisticsGenerator
{
	Params& params;

	struct TaskData
	{
		uint32_t binId;

		TaskData() :
			binId(static_cast<uint32_t>(-1))
		{}
		TaskData(uint32_t binId) :
			binId(binId)
		{}
	};
	std::vector<TaskData> tasksData;
	TasksPool<TaskData> tasksPool;

	std::unique_ptr<kmcdb::MetadataReader> matrixMetadataReader;
	std::unique_ptr<kmcdb::ReaderSortedPlainForListing<uint64_t>> matrixReader;

	using out_kmcdb_value_type = double;

	WritingGatherer<out_kmcdb_value_type> gatherer;
	std::unique_ptr<DumpWriter> normWriter;

	std::unique_ptr<ProgressBar> progress_bar;

	void fillTaskData();


	std::vector<uint8_t> normalizationData;

	const std::vector<int64_t>& correlationPhenotype;
	const std::vector<uint32_t>& differentialAnalysisPhenotype;
	size_t differentialAnalysisNClasses;

	StatisticsToGeneration statisticsToGeneration;

	size_t getMaxNormLineLength() const
	{
		return params.stage1Params.GetKmerLen() + params.mkmcParams.samples.size() * (refresh::numeric_conversion_max_length<out_kmcdb_value_type>() + 1);
	}

	void operator()();

public:
	StatisticsGenerator(Params& params);

	void generateStatisticsParallel();
};

