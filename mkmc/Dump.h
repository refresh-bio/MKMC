#pragma once

#include <vector>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include "parameters.h"
#include "../kmc/kmc_api/kmc_file.h"
#include "KMCFileWrapper.h"
#include "HeapMerge.h"
#include "FileGenerators.h"
#include "TasksPool.h"
#include "Logger.h"
#include "../kmc/kmc_dump/nc_utils.h"
#include "progress_bar.hpp"
#include "Dump.h"
#include "Filter.h"
#include "KmersSamplesStruct.h"

template<unsigned SIZE>
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
	TasksPool<TaskData> tasksPool;

	ProgressBar progress_bar;

	bool allKAreSame(const std::vector<KMCFileWrapper<SIZE>>& samples);
	void openDatabases(std::vector<KMCFileWrapper<SIZE>>& samples, uint32_t binId);
	void fillTaskData();

	uint64_t getNTotInputKmers()
	{
		uint64_t res{};
		for (const auto& x : params.mkmcParams.kmcOutputFiles)
		{
			CKMCFile tmp(true);
			if (!tmp.OpenForListingWithBinOrder(x))
			{
				std::cerr << "Error: cannot open kmc database " << x << "." << std::endl;
				exit(1);
			}
			res += tmp.KmerCount();
		}
		return res;
	}

	template<typename GENERATOR_T, typename Filters>
	void dumpToFile(std::string fileName, uint32_t fileId);

public:
	Dump(const Params& params) :
		params(params), tasksPool(tasksData),
		progress_bar(params.mkmcParams.verbosity_level == 0 ? 0 : getNTotInputKmers(), "Merge", std::cerr, params.mkmcParams.verbosity_level == 0)
	{}

	void dumpToFileParallel();

	void operator()();
};


template<unsigned SIZE>
template<typename GENERATOR_T, typename Filters>
void Dump<SIZE>::dumpToFile(std::string fileName, uint32_t binId)
{
	std::ofstream outputFile(fileName);
	if (!outputFile.is_open())
	{
		std::cerr << "Error: cannot create output file " << params.mkmcParams.outputFilesTemplate << "." << std::endl;
		exit(1);
	}

	std::vector<KMCFileWrapper<SIZE>> samples;
	openDatabases(samples, binId);

	size_t tot_all_kmers{};
	for (const auto& db : samples) {
		tot_all_kmers += db.GetTotKmers();
	}

	ProgressBarUpdater progress_bar_updater(progress_bar, (std::max)(1ull, tot_all_kmers / 100ull));
	uint32_t k = params.stage1Params.GetKmerLen();
	GENERATOR_T fileGenerator(outputFile, params.mkmcParams.samples, params.mkmcParams.count_symbols, k);

	std::vector<size_t> kMersCounts(samples.size());

	Filters filter(params, binId);


	auto do_with_elem_if_exists_init = [&](size_t id, const auto& modifyHeapCallback) -> bool
	{
		modifyHeapCallback(id);
		return !samples[id].Finished();
	};
	auto do_with_elem_if_exists = [&](size_t id, const auto& modifyHeapCallback) -> bool
	{
		assert(!samples[id].Finished());

		samples[id].Next();
		++progress_bar_updater;

		if (samples[id].Finished())
			return false;

		modifyHeapCallback(id);
		return true;
	};

	class HeapComp {
		const std::vector<KMCFileWrapper<SIZE>>& samples;
	public:
		HeapComp(std::vector<KMCFileWrapper<SIZE>>& samples) : samples(samples) {}

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

	CKmer<SIZE> minKmer;

	heap.ProcessElem(do_with_elem_if_exists, [&](size_t elem, size_t id)
		{
			minKmer = samples[elem].First();

			std::fill(kMersCounts.begin(), kMersCounts.end(), 0);
			kMersCounts[id] = samples[elem].FirstCount();
		});

	while (!heap.Empty()) {
		heap.ProcessElem(do_with_elem_if_exists, [&](size_t elem, size_t id)
			{
				const CKmer<SIZE>& curKmer = samples[elem].First();
				if (!(curKmer == minKmer))
				{
					if (filter.keepKMer(KmersSamplesStruct<SIZE>{ minKmer, kMersCounts }))
					{
						fileGenerator.writeKmer(KmersSamplesStruct<SIZE>{ minKmer, kMersCounts });
					}

					minKmer = curKmer;

					std::fill(kMersCounts.begin(), kMersCounts.end(), 0);
				}
				kMersCounts[id] = samples[elem].FirstCount();
			});
	}

	if (filter.keepKMer(KmersSamplesStruct<SIZE>{ minKmer, kMersCounts }))
	{
		fileGenerator.writeKmer(KmersSamplesStruct<SIZE>{ minKmer, kMersCounts });
	}
}


template<unsigned SIZE>
bool Dump<SIZE>::allKAreSame(const std::vector<KMCFileWrapper<SIZE>>& samples)
{
	if (samples.empty())
		return true;
	uint32_t k = samples.front().GetK();
	uint32_t signatureLen = samples.front().GetSignatureLen();
	auto signatureSelectionScheme = samples.front().GetSignatureSelectionScheme();

	for (const auto& sample : samples)
	{
		if (k != sample.GetK())
			return false;
		if (signatureLen != sample.GetSignatureLen())
			return false;
		if (signatureSelectionScheme != sample.GetSignatureSelectionScheme())
			return false;
	}
	return true;
}


template<unsigned SIZE>
void Dump<SIZE>::openDatabases(std::vector<KMCFileWrapper<SIZE>>& samples, uint32_t binId)
{
	for (const std::string& fileName : params.mkmcParams.kmcOutputFiles)
	{
		samples.emplace_back(fileName, binId);
	}

	if (!allKAreSame(samples))
	{
		std::cerr << "Error: KMC databases are not consistent." << std::endl;
		exit(1);
	}
}



template<unsigned SIZE>
inline void Dump<SIZE>::fillTaskData()
{
	tasksData.reserve(params.stage1Params.GetNBins());
	for (uint32_t i = 0; i < params.stage1Params.GetNBins(); ++i)
	{
		tasksData.push_back(TaskData{ i });
	}

	size_t biggestSample = 0;
	uint64_t biggestSampleKmersCount = 0;
	{
		for (size_t i = 0; i < params.mkmcParams.kmcOutputFiles.size(); ++i)
		{
			CKMCFile tmp(true);
			if (!tmp.OpenForListingWithBinOrder(params.mkmcParams.kmcOutputFiles[i]))
			{
				std::cerr << "Error: cannot open kmc database " << params.mkmcParams.kmcOutputFiles[i] << "." << std::endl;
				exit(1);
			}

			if (tmp.KmerCount() > biggestSampleKmersCount)
			{
				biggestSample = i;
				biggestSampleKmersCount = tmp.KmerCount();
			}
		}
	} // close samples

	// sorting is performed in the following manner: first biggest bins are dumped, then smaller; but the sorting is performed basing on the biggest sample only

	std::vector<uint64_t> samplesBeginSize;
	samplesBeginSize.reserve(params.mkmcParams.nKMCBins);
	for (size_t i = 0; i < params.mkmcParams.nKMCBins; ++i)
	{
		KMCFileWrapper<SIZE> currentSample(params.mkmcParams.kmcOutputFiles[biggestSample], i);
		samplesBeginSize.push_back(currentSample.GetTotKmers());
	}

	std::sort(tasksData.begin(), tasksData.end(), [&](const TaskData& a, const TaskData& b) { return samplesBeginSize[a.binId] > samplesBeginSize[b.binId]; });
}



template<unsigned SIZE>
void Dump<SIZE>::dumpToFileParallel()
{
	fillTaskData();

	std::vector<std::thread> threads(params.mkmcParams.nThreads);
	for (uint32_t i_thred = 0; i_thred < params.mkmcParams.nThreads; ++i_thred)
	{
		threads[i_thred] = std::thread([this] { (*this)(); });
	}

	for (std::thread& thread : threads)
	{
		thread.join();
	}
}


template<unsigned SIZE>
void Dump<SIZE>::operator()()
{
	TaskData taskData;
	using ParameterizedKmersSamplesStruct = KmersSamplesStruct<SIZE>;

	if (params.mkmcParams.outputFileType == OutputFileType::Matrix)
		while (tasksPool.getTask(taskData))
		{
			if (params.filterParams.filterKmersSequences)
			{
				using Filters = FilterCountThreshold<ParameterizedKmersSamplesStruct, FilterSequences<ParameterizedKmersSamplesStruct>>;
				dumpToFile<MatrixFileGenerator, Filters>(params.mkmcParams.outputFiles[taskData.binId], taskData.binId);
			}
			else
			{
				using Filters = FilterCountThreshold<ParameterizedKmersSamplesStruct>;
				dumpToFile<MatrixFileGenerator, Filters>(params.mkmcParams.outputFiles[taskData.binId], taskData.binId);
			}
		}
	else if (params.mkmcParams.outputFileType == OutputFileType::FASTA)
		while (tasksPool.getTask(taskData))
		{
			if (params.filterParams.filterKmersSequences)
			{
				using Filters = FilterCountThreshold<ParameterizedKmersSamplesStruct, FilterSequences<ParameterizedKmersSamplesStruct>>;
				dumpToFile<FASTAFileGenerator, Filters>(params.mkmcParams.outputFiles[taskData.binId], taskData.binId);
			}
			else
			{
				using Filters = FilterCountThreshold<ParameterizedKmersSamplesStruct>;
				dumpToFile<FASTAFileGenerator, Filters>(params.mkmcParams.outputFiles[taskData.binId], taskData.binId);
			}
		}
	else
		assert(false);
}
