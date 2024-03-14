#pragma once

#include <vector>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include "parameters.h"
#include "../kmc/kmc_api/kmc_file.h"
#include "KMCFileWrapper.h"
#include "HeapMerge.h"
#include "FileGenerators.h"
#include "Logger.h"



class Dump
{
	const Params& params;

	bool allKAreSame(const std::vector<KMCFileWrapper>& samples);
	void openDatabases(std::vector<KMCFileWrapper>& samples);

public:
	Dump(const Params& params) :
		params(params)
	{}

	void dumpToStd();
	template<typename GENRATOR_T>
	void dumpToFile();
};



template<typename GENRATOR_T>
void Dump::dumpToFile()
{
	std::ofstream outputFile(params.mkmcParams.outputFile);
	if (!outputFile.is_open())
	{
		std::cerr << "Error: cannot create output file " << params.mkmcParams.outputFile << "." << std::endl;
		exit(1);
	}

	std::vector<KMCFileWrapper> samples;
	openDatabases(samples);

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

	BinaryHeapMergeStreams<size_t, HeapComp> heap(samples.size(), do_with_elem_if_exists_init, HeapComp(samples));

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
