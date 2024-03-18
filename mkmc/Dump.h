#pragma once

#include <vector>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <algorithm>
#include "parameters.h"
#include "../kmc/kmc_api/kmc_file.h"
#include "KMCFileWrapper.h"
#include "HeapMerge.h"
#include "FileGenerators.h"
#include "TasksPool.h"
#include "Logger.h"



class Dump
{
	const Params& params;

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

	bool allKAreSame(const std::vector<KMCFileWrapper>& samples);
	void openDatabases(std::vector<KMCFileWrapper>& samples, uint32_t binId);

	template<typename GENRATOR_T>
	void dumpToFile(std::string fileName, uint32_t fileId);

public:
	Dump(const Params& params) :
		params(params), tasksPool(tasksData)
	{}

	template<typename GENRATOR_T>
	void dumpToFileParallel();

	TasksPool<TaskData> tasksPool;

	template<typename GENRATOR_T>
	void operator()();
};



template<typename GENRATOR_T>
void Dump::dumpToFile(std::string fileName, uint32_t binId)
{
	std::ofstream outputFile(fileName);
	if (!outputFile.is_open())
	{
		std::cerr << "Error: cannot create output file " << params.mkmcParams.outputFilesTemplate << "." << std::endl;
		exit(1);
	}

	std::vector<KMCFileWrapper> samples;
	openDatabases(samples, binId);

	size_t tot_all_kmers{};
	for (const auto& db : samples) {
		tot_all_kmers += db.GetTotKmers();
	}
	PercentProgress progress(tot_all_kmers, params.mkmcParams.verbosity_level > 0);

	uint32_t k = params.stage1Params.GetKmerLen();
	GENRATOR_T fileGenerator(outputFile, params.mkmcParams.samples, params.mkmcParams.count_symbols, k);

	std::vector<size_t> kMersCounts(samples.size());
	Filter filter(params);


	auto do_with_elem_if_exists_init = [&](size_t id, const auto& modifyHeapCallback) -> bool
	{
		modifyHeapCallback(id);
		return !samples[id].Finished();
	};
	auto do_with_elem_if_exists = [&](size_t id, const auto& modifyHeapCallback) -> bool
	{
		assert(!samples[id].Finished());

		samples[id].Next();
		progress.NotifyProgress(1);

		if (samples[id].Finished())
			return false;

		modifyHeapCallback(id);
		return true;
	};

	class HeapComp {
		const std::vector<KMCFileWrapper>& samples;
	public:
		HeapComp(std::vector<KMCFileWrapper>& samples) : samples(samples) {}

		bool operator()(const size_t a, const size_t b) const
		{
			return samples[b].First() < samples[a].First();
		}
	};

	std::vector<size_t> streams_to_merge;
	for (size_t i = 0; i < samples.size(); ++i)
	{
		if (!samples[i].Finished())
			streams_to_merge.push_back(i);
	}

	BinaryHeapMergeStreams<size_t, HeapComp> heap(streams_to_merge, do_with_elem_if_exists_init, HeapComp(samples));

	if (heap.Empty())
		return;

	KMCFileWrapper::kmer_t minKmer;

	heap.ProcessElem(do_with_elem_if_exists, [&](size_t elem, size_t id)
		{
			minKmer = samples[elem].First();

			std::fill(kMersCounts.begin(), kMersCounts.end(), 0);
			kMersCounts[id] = samples[elem].FirstCount();
		});

	while (!heap.Empty()) {
		heap.ProcessElem(do_with_elem_if_exists, [&](size_t elem, size_t id)
			{
				const KMCFileWrapper::kmer_t& curKmer = samples[elem].First();
				if (!(curKmer == minKmer))
				{
					if (filter.keepKMer(kMersCounts))
					{
						fileGenerator.writeKmer(minKmer, kMersCounts);
					}

					minKmer = curKmer;

					std::fill(kMersCounts.begin(), kMersCounts.end(), 0);
				}
				kMersCounts[id] = samples[elem].FirstCount();
			});
	}

	if (filter.keepKMer(kMersCounts))
	{
		fileGenerator.writeKmer(minKmer, kMersCounts);
	}
}



template<typename GENRATOR_T>
void Dump::dumpToFileParallel()
{
	tasksData.reserve(params.stage1Params.GetNBins());
	for (uint32_t i = 0; i < params.stage1Params.GetNBins(); ++i)
	{
		tasksData.push_back(TaskData{ i });
	}

	std::vector<std::thread> threads(params.stage1Params.GetNThreads());
	for (uint32_t i_thred = 0; i_thred < params.mkmcParams.nKMCWorkers; ++i_thred)
	{
		threads[i_thred] = std::thread([this] { (*this).operator()<GENRATOR_T>(); });
	}

	for (std::thread& thread : threads)
	{
		thread.join();
	}
}

template<typename GENRATOR_T>
void Dump::operator()()
{
	TaskData taskData;
	while (tasksPool.getTask(taskData))
	{
		dumpToFile<GENRATOR_T>(params.mkmcParams.outputFiles[taskData.binId], taskData.binId);
	}
}
