#pragma once

#include "FileGenerators.h"
#include <vector>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <mutex>
#include "parameters.h"
#include "kmc_dump/nc_utils.h"
#include "KMCFileWrapper.h"
#include "HeapMerge.h"
#include "TasksPool.h"
#include "Logger.h"
#define NOMINMAX
#include "progress_bar.hpp"
#include "Filter.h"
#include "MatrixStats.h"
#include "KmersSamplesStruct.h"
#include "kmcdb/kmcdb.h"
#include "refresh/statistics/lib/statistics_normalization.h"

template<unsigned SIZE>
class Dump
{
	const Params& params;

	StatisticsParams::NormalizationLearning normalizationLearning;
	std::mutex normalizationLearningMutex;

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

	std::vector<std::unique_ptr<kmcdb::MetadataReader>> samplesMetadata;
	std::vector<std::unique_ptr<kmcdb::ReaderSortedWithLUTForListing<uint64_t>>> samplesReaders;

	std::unique_ptr<kmcdb::MetadataReader> sequencesToFilterMetadataReader;
	std::unique_ptr<kmcdb::ReaderSortedWithLUTForListing<uint64_t>> sequencesToFilterReader;

	std::unique_ptr<ProgressBar> progress_bar;

	bool inputIsConsistent();
	void fillTaskData();
	void serializeNormalizationAndDump();

	template<typename Generators_T, typename Filters_T>
	void dumpToFile(uint32_t binId, StatisticsParams::NormalizationLearning& normalizationLearning, Generators_T& fileGenerators, kmcdb::BinReaderSortedWithLUTForListing<uint64_t>* bin);

public:
	Dump(const Params& params) :
		params(params), tasksPool(tasksData)
	{
		normalizationLearning.register_method(StatisticsParams::NormalizationMethod::frequency_count);
		normalizationLearning.register_method(StatisticsParams::NormalizationMethod::quantile);
		normalizationLearning.set_no_series(params.mkmcParams.samples.size());
		normalizationLearning.initialize();
	}

	void dumpToFileParallel();

	template<typename Generators_T>
	void operator()();
};



template<unsigned SIZE>
template<typename Generators_T, typename Filters_T>
void Dump<SIZE>::dumpToFile(uint32_t binId, StatisticsParams::NormalizationLearning& normalizationLearning, Generators_T& fileGenerators, kmcdb::BinReaderSortedWithLUTForListing<uint64_t>* bin)
{
	std::vector<KMCFileWrapper<SIZE>> samples;
	for (size_t sample_id = 0; sample_id < params.mkmcParams.kmcOutputFiles.size(); ++sample_id)
		samples.emplace_back(samplesReaders[sample_id]->GetBin(binId));

	size_t tot_all_kmers{};
	for (const auto& db : samples) {
		tot_all_kmers += db.GetTotKmers();
	}

	ProgressBarUpdater progress_bar_updater(*progress_bar, (std::max)(1ull, tot_all_kmers / 100ull));

	std::vector<uint64_t> kMersCounts(samples.size());

	Filters_T filter(params, bin);

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

	kmcdb::CKmer<SIZE> minKmer;

	heap.ProcessElem(do_with_elem_if_exists, [&](size_t elem, size_t id)
		{
			minKmer = samples[elem].First();

			std::fill(kMersCounts.begin(), kMersCounts.end(), 0);
			kMersCounts[id] = samples[elem].FirstCount();
		});

	while (!heap.Empty()) {
		heap.ProcessElem(do_with_elem_if_exists, [&](size_t elem, size_t id)
			{
				const kmcdb::CKmer<SIZE>& curKmer = samples[elem].First();
				if (!(curKmer == minKmer))
				{
					if (filter.keepKMer(KmersSamplesStruct<SIZE>{ minKmer, kMersCounts }))
					{
						fileGenerators.writeKmer(KmersSamplesStruct<SIZE>{ minKmer, kMersCounts });
						normalizationLearning.add_entry(kMersCounts);
					}

					minKmer = curKmer;

					std::fill(kMersCounts.begin(), kMersCounts.end(), 0);
				}
				kMersCounts[id] = samples[elem].FirstCount();
			});
	}

	if (filter.keepKMer(KmersSamplesStruct<SIZE>{ minKmer, kMersCounts }))
	{
		fileGenerators.writeKmer(KmersSamplesStruct<SIZE>{ minKmer, kMersCounts });
		normalizationLearning.add_entry(kMersCounts);
	}
}


template<unsigned SIZE>
bool Dump<SIZE>::inputIsConsistent()
{
	if (samplesMetadata.empty())
		return true;
	uint32_t k = samplesMetadata.front()->GetConfig().kmer_len;
	uint64_t signatureLen = samplesMetadata.front()->GetConfig().signature_len;
	auto signatureSelectionScheme = samplesMetadata.front()->GetConfig().signature_selection_scheme;
	auto signatureToBinMapping = samplesMetadata.front()->GetConfig().signature_to_bin_mapping;
	auto num_bins = samplesMetadata.front()->GetConfig().num_bins;

	for (const auto& sample : samplesMetadata)
	{
		if (k != sample->GetConfig().kmer_len)
			return false;
		if (signatureLen != sample->GetConfig().signature_len)
			return false;
		if (signatureSelectionScheme != sample->GetConfig().signature_selection_scheme)
			return false;
		if (signatureToBinMapping != sample->GetConfig().signature_to_bin_mapping)
			return false;
		if (num_bins != sample->GetConfig().num_bins)
			return false;
	}
	return true;
}

template<unsigned SIZE>
inline void Dump<SIZE>::fillTaskData()
{
	tasksData.reserve(params.stage1Params.GetNBins());
	for (uint32_t i = 0; i < params.stage1Params.GetNBins(); ++i)
	{
		tasksData.push_back(TaskData{ i });
	}

	uint64_t biggestSampleKmersCount = 0;

	std::vector<uint64_t> samplesBeginSize(params.mkmcParams.nKMCBins);

	samplesMetadata.reserve(params.stage1Params.GetNBins());
	samplesReaders.reserve(params.stage1Params.GetNBins());

	uint64_t totKmersAllSamples = 0;

	for (size_t i = 0; i < params.mkmcParams.kmcOutputFiles.size(); ++i)
	{
		try
		{
			samplesMetadata.emplace_back(std::make_unique<kmcdb::MetadataReader>(params.mkmcParams.kmcOutputFiles[i], true));
			kmcdb::MetadataReader& metadata_reader = *samplesMetadata.back();
			samplesReaders.emplace_back(std::make_unique<kmcdb::ReaderSortedWithLUTForListing<uint64_t>>(metadata_reader));
			kmcdb::ReaderSortedWithLUTForListing<uint64_t>& reader = *samplesReaders.back();
			uint64_t totKmers = 0;
			for (uint32_t bin_id = 0; bin_id < metadata_reader.GetConfig().num_bins; ++bin_id)
				totKmers += reader.GetBin(bin_id)->GetBinMetadata().total_kmers;

			totKmersAllSamples += totKmers;
			if (totKmers > biggestSampleKmersCount)
			{
				biggestSampleKmersCount = totKmers;

				for (uint32_t bin_id = 0; bin_id < metadata_reader.GetConfig().num_bins; ++bin_id)
					samplesBeginSize[bin_id] = reader.GetBin(bin_id)->GetBinMetadata().total_kmers;
			}
		}
		catch (const std::runtime_error& ex)
		{
			std::cerr << "Error: " << ex.what() << std::endl;
			exit(1);
		}
	}

	if (!inputIsConsistent())
	{
		std::cerr << "Error: KMC databases are not consistent." << std::endl;
		exit(1);
	}

	progress_bar = std::make_unique<ProgressBar>(params.mkmcParams.verbosity_level == 0 ? 0 : totKmersAllSamples, "Dumping", std::cerr, params.mkmcParams.verbosity_level == 0);

	// sorting is performed in the following manner: first biggest bins are dumped, then smaller; but the sorting is performed basing on the biggest sample only
	std::sort(tasksData.begin(), tasksData.end(), [&](const TaskData& a, const TaskData& b) { return samplesBeginSize[a.binId] > samplesBeginSize[b.binId]; });
}

template<unsigned SIZE>
void Dump<SIZE>::serializeNormalizationAndDump()
{
	std::vector<uint8_t> frequencyNormalizationData, quantileNormalizationData;

	normalizationLearning.serialize(StatisticsParams::NormalizationMethod::frequency_count, frequencyNormalizationData);
	normalizationLearning.serialize(StatisticsParams::NormalizationMethod::quantile, quantileNormalizationData);

	MatrixStatsWriter stats_writer(params.mkmcParams.normStatsBinFile);
	stats_writer.Add(params.statisticsParams.normFrequencyStreamName, frequencyNormalizationData);
	stats_writer.Add(params.statisticsParams.normQuantileStreamName, quantileNormalizationData);
}

template<unsigned SIZE>
void Dump<SIZE>::dumpToFileParallel()
{
	fillTaskData();

	if (params.filterParams.filterKmersSequences)
	{
		sequencesToFilterMetadataReader = std::make_unique<kmcdb::MetadataReader>(params.filterParams.kmersSequencesToFilterOutDB, true);
		sequencesToFilterReader = std::make_unique<kmcdb::ReaderSortedWithLUTForListing<uint64_t>>(*sequencesToFilterMetadataReader);
	}

	kmcdb::Config config;
	{
		config.num_bins = samplesMetadata.front()->GetConfig().num_bins;
		config.signature_len = samplesMetadata.front()->GetConfig().signature_len;
		config.signature_selection_scheme = samplesMetadata.front()->GetConfig().signature_selection_scheme;
		config.signature_to_bin_mapping = samplesMetadata.front()->GetConfig().signature_to_bin_mapping;
		config.kmer_len = samplesMetadata.front()->GetConfig().kmer_len;
		config.num_samples = samplesMetadata.size();
		config.num_bytes_single_value = samplesMetadata.front()->GetConfig().num_bytes_single_value;
		for (size_t i = 1; i < samplesMetadata.size(); ++i)
			if (samplesMetadata[i]->GetConfig().num_bytes_single_value > config.num_bytes_single_value)
				config.num_bytes_single_value = samplesMetadata[i]->GetConfig().num_bytes_single_value;
	}

	std::vector<std::thread> threads(params.mkmcParams.nThreads);

	if (params.mkmcParams.outputFileTypes.empty())
	{
		using Generators_T = PerformGenerate<BinFileGenerator>;
		Generators_T::initWriters(params, config);
		for (uint32_t i_thred = 0; i_thred < params.mkmcParams.nThreads; ++i_thred)
		{
			threads[i_thred] = std::thread([this] { this->operator() < Generators_T > (); });
		}

		for (std::thread& thread : threads)
		{
			thread.join();
		}
		Generators_T::closeWriters();
	}
	else if (params.mkmcParams.outputFileTypes.size() == 2)
	{
		using Generators_T = PerformGenerate<BinFileGenerator, MatrixFileGenerator, FASTAFileGenerator>;
		Generators_T::initWriters(params, config);
		for (uint32_t i_thred = 0; i_thred < params.mkmcParams.nThreads; ++i_thred)
		{
			threads[i_thred] = std::thread([this] { this->operator() < Generators_T > (); });
		}

		for (std::thread& thread : threads)
		{
			thread.join();
		}
		Generators_T::closeWriters();
	}
	else if (params.mkmcParams.outputFileTypes.front() == OutputFileType::Matrix)
	{
		using Generators_T = PerformGenerate<BinFileGenerator, MatrixFileGenerator>;
		Generators_T::initWriters(params, config);
		for (uint32_t i_thred = 0; i_thred < params.mkmcParams.nThreads; ++i_thred)
		{
			threads[i_thred] = std::thread([this] { this->operator() < Generators_T > (); });
		}

		for (std::thread& thread : threads)
		{
			thread.join();
		}
		Generators_T::closeWriters();
	}
	else
	{
		using Generators_T = PerformGenerate<BinFileGenerator, FASTAFileGenerator>;
		Generators_T::initWriters(params, config);
		for (uint32_t i_thred = 0; i_thred < params.mkmcParams.nThreads; ++i_thred)
		{
			threads[i_thred] = std::thread([this] { this->operator() < Generators_T > (); });
		}
		
		for (std::thread& thread : threads)
		{
			thread.join();
		}
		Generators_T::closeWriters();
	}

	if (params.statisticsParams.generateNormalization || params.statisticsParams.generateEntropy || !params.statisticsParams.classificationMethods.empty())
		serializeNormalizationAndDump();
}


template<unsigned SIZE>
template<typename Generators_T>
void Dump<SIZE>::operator()()
{
	Generators_T fileGenerators(params);
	TaskData taskData;
	using ParameterizedKmersSamplesStruct = KmersSamplesStruct<SIZE>;

	while (tasksPool.getTask(taskData))
	{
		StatisticsParams::NormalizationLearning currentBinNormalizationLearnings;
		currentBinNormalizationLearnings.register_method(StatisticsParams::NormalizationMethod::frequency_count);
		currentBinNormalizationLearnings.register_method(StatisticsParams::NormalizationMethod::quantile);
		currentBinNormalizationLearnings.set_no_series(params.mkmcParams.samples.size());
		currentBinNormalizationLearnings.initialize();

		fileGenerators.setBinId(taskData.binId);

		if (params.filterParams.filterKmersSequences)
		{
			assert(sequencesToFilterReader);
			using Filters = PerformFilter<FilterCountThreshold<ParameterizedKmersSamplesStruct>, FilterSequences<ParameterizedKmersSamplesStruct>>;
			dumpToFile<Generators_T, Filters>(taskData.binId, currentBinNormalizationLearnings, fileGenerators, sequencesToFilterReader->GetBin(taskData.binId));
		}
		else
		{
			using Filters = PerformFilter<FilterCountThreshold<ParameterizedKmersSamplesStruct>>;
			dumpToFile<Generators_T, Filters>(taskData.binId, currentBinNormalizationLearnings, fileGenerators, nullptr);
		}
		normalizationLearningMutex.lock();
		normalizationLearning.merge_with(&currentBinNormalizationLearnings, &currentBinNormalizationLearnings + 1);
		normalizationLearningMutex.unlock();
	}
}
