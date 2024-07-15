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

	struct
	{
		bool pearson = false;
		bool spearman = false;
		bool kendall = false;

		bool entropy = false;

		bool differentialAnalysis = false; // logical sum of the following ones

		bool tTest = false;
		bool snr = false;
		bool wilcoxonRankSum = false;

		bool dids = false;
		bool anova = false;
	} statisticsToGeneration;

	struct
	{
		std::unique_ptr<DumpWriter> pearson;
		std::unique_ptr<DumpWriter> spearman;
		std::unique_ptr<DumpWriter> kendall;

		std::unique_ptr<DumpWriter> entropy;

		std::unique_ptr<DumpWriter> tTest;
		std::unique_ptr<DumpWriter> snr;
		std::unique_ptr<DumpWriter> wilcoxonRankSum;

		std::unique_ptr<DumpWriter> dids;
		std::unique_ptr<DumpWriter> anova;
	} writers;

public:
	StatisticsGenerator(Params& params);

	void generateStatisticsParallel();
};

