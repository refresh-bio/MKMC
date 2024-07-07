#pragma once

#include <fstream>
#include <vector>
#include <string>
#include "parameters.h"
#include "TasksPool.h"
#include "refresh/statistics/lib/statistics.h"
#define NOMINMAX
#include "progress_bar.hpp"
#include "kmcdb/kmcdb.h"


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

	std::unique_ptr<kmcdb::WriterSortedPlain<out_kmcdb_value_type>> kmcdbWriter;

	std::unique_ptr<ProgressBar> progress_bar;

	void fillTaskData();


	std::vector<uint8_t> normalizationData;

	const std::vector<int64_t>& correlationPhenotype;
	const std::vector<uint32_t>& differentialAnalysisPhenotype;
	size_t differentialAnalysisNClasses;

	void operator()();

public:
	StatisticsGenerator(Params& params) :
		params(params),
		tasksPool(tasksData),
		correlationPhenotype(params.phenotypes.correlationPhenotype.getPhenotype()),
		differentialAnalysisPhenotype(params.phenotypes.differentialAnalysisPhenotype.getMappedPhenotype()),
		differentialAnalysisNClasses(params.phenotypes.differentialAnalysisPhenotype.getClassesNumber())
	{}

	void generateStatisticsParallel();
};

