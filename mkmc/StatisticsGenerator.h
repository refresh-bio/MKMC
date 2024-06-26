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

	uint64_t totAllKmers;
	ProgressBar progress_bar;

	void fillTaskData();

	void readPhenotype(std::vector<int>& phenotype);
	template<typename T>
	void readDump(std::vector<T>& normalizationData, std::string normalizationFileName);

	bool getLine(std::ifstream& stream, std::string& kmerSequence, std::vector<uint64_t>& counts)
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

	uint64_t getNTotInputKmers()
	{
		std::vector<uint64_t> nOutputKmersPerBin;
		readDump(nOutputKmersPerBin, params.statisticsParams.statsNOutputKmers);
		uint64_t nKmers = 0;
		for (auto a : nOutputKmersPerBin)
			nKmers += a;
		return nKmers;
	}

	std::vector<uint8_t> normalizationData;
	std::vector<int> phenotype;

	void operator()();

public:
	StatisticsGenerator(Params& params) :
		params(params),
		tasksPool(tasksData),
		totAllKmers(getNTotInputKmers()),
		progress_bar(params.mkmcParams.verbosity_level == 0 ? 0 : totAllKmers, "Computing statistics", std::cerr, params.mkmcParams.verbosity_level == 0)
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
