#pragma once

#include <fstream>
#include <vector>
#include <string>
#include "parameters.h"
#include "TasksPool.h"
#include "lib/statistics/lib/statistics.h"



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

	void readPhenotype(std::vector<int>& phenotype);
	void readNormalizationDump(std::vector<uint8_t>& normalizationData, std::string normalizationFileName);

	bool getLine(std::ifstream& stream, std::string& kmerSequence, std::vector<size_t>& counts)
	{
		stream >> kmerSequence;
		if (stream.eof())
			return false;
		for (auto& count : counts)
		{
			stream >> count;
		}
		return true;
	}
	void putLine(std::ofstream& stream, std::string& kmerSequence, const std::vector<double>& counts)
	{
		stream << kmerSequence << '\t';
		for (auto count : counts)
		{
			stream << count << '\t';
		}
		stream << '\n';
	}

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
