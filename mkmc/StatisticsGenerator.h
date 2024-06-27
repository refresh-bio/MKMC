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

	void readPhenotype(std::vector<int>& phenotype);
	template<typename T>
	void readDump(std::vector<T>& normalizationData, std::string normalizationFileName);

	std::vector<uint8_t> normalizationData;
	std::vector<int> phenotype;

	void operator()();

public:
	StatisticsGenerator(Params& params) :
		params(params),
		tasksPool(tasksData)
	{}

	void generateStatisticsParallel();
};

template<typename T>
void StatisticsGenerator::readDump(std::vector<T>& data, std::string fileName)
{
	std::ifstream file(fileName, std::ios::binary);
	if (!file.is_open())
	{
		std::cerr << "Error: cannot open " << fileName << "." << std::endl;
		exit(1);
	}
	size_t nElements;
	file.read(reinterpret_cast<char*>(&nElements), sizeof(size_t));
	data.resize(nElements);
	file.read(reinterpret_cast<char*>(data.data()), nElements * sizeof(T));
}
