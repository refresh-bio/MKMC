#pragma once

#include <vector>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include "parameters.h"
#include "../kmc/kmc_api/kmc_file.h"
#include "KMCFileWrapper.h"
#include "FileGenerators.h"



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



//#define POP_HEAP
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

	uint32_t k = params.stage1Params.GetKmerLen();
	GENRATOR_T fileGenerator(outputFile, params.mkmcParams.inputFiles, params.mkmcParams.count_symbols, k);

	std::vector<size_t> kMersCounts(samples.size());
	Filter filter(params);

	class HeapComp {
		std::vector<KMCFileWrapper>& samples;
	public:
		HeapComp(std::vector<KMCFileWrapper>& samples) : samples(samples) {}

		bool operator()(size_t a, size_t b)
		{
			// aFinished && !bFinished -> true
			// bFinished && !aFinished -> false
			// aFinished && bFinished -> false

			if (samples[a].Finished())
			{
				return !samples[b].Finished();
			}
			if (samples[b].Finished())
			{
				return false;
			}
			return samples[b].First() < samples[a].First();
		}
	};
	std::vector<size_t> kmersHeap(samples.size());
	std::iota(kmersHeap.begin(), kmersHeap.end(), 0);
	std::make_heap(kmersHeap.begin(), kmersHeap.end(), HeapComp(samples));

	while (true)
	{
		size_t minId = kmersHeap.front();
#ifndef POP_HEAP
		if (samples[minId].Finished())
			break;
#endif

		std::fill(kMersCounts.begin(), kMersCounts.end(), 0);

		auto minKmer = samples[minId].First();

		kMersCounts[minId] = samples[minId].FirstCount();
		samples[minId].Next();

		while (true)
		{
			std::pop_heap(kmersHeap.begin(), kmersHeap.end(), HeapComp(samples));

#ifdef POP_HEAP
			if (samples[kmersHeap.back()].Finished())
			{
				kmersHeap.pop_back();
				if (kmersHeap.empty())
					break;
	}
			else
				std::push_heap(kmersHeap.begin(), kmersHeap.end(), HeapComp(samples));
#else
			std::push_heap(kmersHeap.begin(), kmersHeap.end(), HeapComp(samples));
#endif

			size_t curId = kmersHeap.front();
			KMCFileWrapper& curKmer = samples[curId];
			if (curKmer.Finished() || !(curKmer.First() == minKmer))
			{
				break;
			}

			kMersCounts[curId] = curKmer.FirstCount();
			samples[curId].Next();
}

		if (filter.keepKMer(kMersCounts))
		{
			fileGenerator.writeKmer(minKmer, kMersCounts);
		}
	}
}