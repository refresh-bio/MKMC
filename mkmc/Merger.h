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

inline void StoreTotCnt(std::vector<std::vector<uint64>>& tot_cnts,
	const std::vector<std::string>& sample_names,
	const Params& params)
{
	//reduce
	for (size_t i = 1; i < tot_cnts.size(); ++i)
		for (size_t sample_id = 0; sample_id < tot_cnts.front().size(); ++sample_id)
			tot_cnts[0][sample_id] += tot_cnts[i][sample_id];

	TextFileWriter writer(params.mkmcParams.outputFileTotCnt, false);
	writer.StoreHeader(sample_names, "measure");

	std::string first_col = "tot_cnt";

	MatrixOutputBuffer<uint64_t> out(writer, first_col.length(), tot_cnts.front().size());
	out.StoreKmer(first_col, tot_cnts[0]);
}

template<unsigned SIZE>
class Merger
{
	const Params& params;

	StatisticsParams::NormalizationLearning normalizationLearning;
	std::mutex normalizationLearningMutex;

	bool matrixNotEmpty = false;

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

	void openReadersAndVerifySamples(/* out */uint64_t& totKmersAllSamples, /* out */size_t& biggestSample);
	void fillTaskData(const size_t biggestSample);
	void serializeNormalizationAndSave();

	template<typename PerformGenerate_T>
	void callThreads(const kmcdb::Config& config, std::vector<std::vector<uint64_t>>& tot_cnts, uint64_t totKmersAllSamples);

	template<unsigned I = static_cast<unsigned>(OutputFileType::_N), typename... Generators_T>
	class GeneratorVecToTemplate {
	public:
		template<typename Callback_T>
		static void callTemplateFunction(const std::vector<OutputFileType>& generators, const Callback_T& callback) {
			if (I > generators.size()) {
				GeneratorVecToTemplate<I - 1, Generators_T...>::callTemplateFunction(generators, callback);
				return;
			}

			switch (generators[generators.size() - I]) {
			case OutputFileType::FASTA:
				GeneratorVecToTemplate<I - 1, Generators_T..., FASTAFileGenerator>::callTemplateFunction(generators, callback); break;
			case OutputFileType::Matrix:
				GeneratorVecToTemplate<I - 1, Generators_T..., MatrixFileGenerator>::callTemplateFunction(generators, callback); break;
			default:
				assert(false);
			}
		}
	};

	template<typename... Generators_T>
	class GeneratorVecToTemplate<0, Generators_T...> {
	public:
		template<typename Callback_T>
		static void callTemplateFunction(const std::vector<OutputFileType>& generators, const Callback_T& callback) {
			callback.template operator() < Generators_T... > ();
		}
	};

	template<typename PerformGenerate_T, typename Filters_T>
	bool mergeToGenerators(uint32_t binId, std::vector<uint64_t>& tot_cnts, StatisticsParams::NormalizationLearning& currentBinNormalizationLearning, PerformGenerate_T& fileGenerators, kmcdb::BinReaderSortedWithLUTForListing<uint64_t>* bin);

public:
	Merger(const Params& params) :
		params(params), tasksPool(tasksData)
	{
		if (params.statisticsParams.normalizationMethod == StatisticsParams::NormalizationMethod::deseq2 || params.statisticsParams.learnDeseq2)
			normalizationLearning.register_method(StatisticsParams::NormalizationMethod::deseq2);
		normalizationLearning.register_method(StatisticsParams::NormalizationMethod::frequency_count);
		normalizationLearning.register_method(StatisticsParams::NormalizationMethod::quantile);
		normalizationLearning.set_no_series(params.mkmcParams.samples.size());
		normalizationLearning.initialize();
	}

	bool mergeParallel();

	template<typename PerformGenerate_T>
	void operator()(std::vector<uint64_t>& tot_cnts);
};



template<unsigned SIZE>
template<typename PerformGenerate_T, typename Filters_T>
bool Merger<SIZE>::mergeToGenerators(uint32_t binId, std::vector<uint64_t>& tot_cnts, StatisticsParams::NormalizationLearning& currentBinNormalizationLearning, PerformGenerate_T& fileGenerators, kmcdb::BinReaderSortedWithLUTForListing<uint64_t>* bin)
{
	std::vector<KMCFileWrapper<SIZE>> samples;
	for (size_t sample_id = 0; sample_id < params.mkmcParams.kmcOutputFiles.size(); ++sample_id)
		samples.emplace_back(samplesReaders[sample_id]->GetBin(binId));

	size_t totAllBinKmers{};
	for (const auto& db : samples) {
		totAllBinKmers += db.GetTotKmers();
	}

	ProgressBarUpdater progress_bar_updater(*progress_bar, (std::max)(1ull, totAllBinKmers / 100ull));

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
		return false;

	uint64_t outputKmerId = 0;
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
						for (size_t i = 0; i < kMersCounts.size(); ++i)
							tot_cnts[i] += kMersCounts[i];

						fileGenerators.writeKmer(KmersSamplesStruct<SIZE>{ minKmer, kMersCounts }, outputKmerId);
						++outputKmerId;
						currentBinNormalizationLearning.add_entry(kMersCounts);
					}

					minKmer = curKmer;

					std::fill(kMersCounts.begin(), kMersCounts.end(), 0);
				}
				kMersCounts[id] = samples[elem].FirstCount();
			});
	}

	if (filter.keepKMer(KmersSamplesStruct<SIZE>{ minKmer, kMersCounts }))
	{
		for (size_t i = 0; i < kMersCounts.size(); ++i)
			tot_cnts[i] += kMersCounts[i];

		fileGenerators.writeKmer(KmersSamplesStruct<SIZE>{ minKmer, kMersCounts }, outputKmerId);
		++outputKmerId;
		currentBinNormalizationLearning.add_entry(kMersCounts);
	}
	return outputKmerId != 0;
}


template<unsigned SIZE>
bool Merger<SIZE>::inputIsConsistent()
{
	if (samplesMetadata.empty())
		return true;
	uint64_t k = samplesMetadata.front()->GetConfig().kmer_len;
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
void Merger<SIZE>::openReadersAndVerifySamples(/* out */uint64_t& totKmersAllSamples, /* out */size_t& biggestSample)
{
	samplesMetadata.reserve(params.mkmcParams.kmcOutputFiles.size());
	samplesReaders.reserve(params.mkmcParams.kmcOutputFiles.size());

	uint64_t biggestSampleKmersCount = 0;

	std::vector<size_t> emptySamplesIndices;

	for (size_t i = 0; i < params.mkmcParams.kmcOutputFiles.size(); ++i) // iterate on samples
	{
		try
		{
			samplesMetadata.emplace_back(std::make_unique<kmcdb::MetadataReader>(params.mkmcParams.kmcOutputFiles[i], true));
			kmcdb::MetadataReader& metadata_reader = *samplesMetadata.back();
			samplesReaders.emplace_back(std::make_unique<kmcdb::ReaderSortedWithLUTForListing<uint64_t>>(metadata_reader));
			kmcdb::ReaderSortedWithLUTForListing<uint64_t>& reader = *samplesReaders.back();

			uint64_t totSampleKmers = 0;
			for (uint32_t bin_id = 0; bin_id < metadata_reader.GetConfig().num_bins; ++bin_id)
				totSampleKmers += reader.GetBin(bin_id)->GetBinMetadata().total_kmers;

			if (totSampleKmers == 0)
				emptySamplesIndices.push_back(i);

			totKmersAllSamples += totSampleKmers;
			if (totSampleKmers > biggestSampleKmersCount)
			{
				biggestSampleKmersCount = totSampleKmers;
				biggestSample = i;
			}
		}
		catch (const std::runtime_error& ex)
		{
			Logger::Inst().Log(std::string("Error: ") + ex.what());
			exit(1);
		}
	}

	if (totKmersAllSamples == 0)
	{
		Logger::Inst().Log("Error: no k-mers present in samples; input files are empty or --ci and --cx parameters are too strict.");
		exit(1);
	}
	else
		for (auto i : emptySamplesIndices)
			Logger::Inst().Log("Warning: a sample " + params.mkmcParams.samples[i].name + " has no k-mers; its input files are empty or --ci and --cx parameters are too strict.");

	if (!inputIsConsistent())
	{
		Logger::Inst().Log("Error: KMC databases are not consistent. Please contact the authors.");
		exit(1);
	}
}

template<unsigned SIZE>
void Merger<SIZE>::fillTaskData(const size_t biggestSample)
{
	const uint64_t numBins = samplesMetadata.front()->GetConfig().num_bins;

	std::vector<uint64_t> samplesBeginSize(numBins);
	for (uint32_t bin_id = 0; bin_id < numBins; ++bin_id)
		samplesBeginSize[bin_id] = samplesReaders[biggestSample]->GetBin(bin_id)->GetBinMetadata().total_kmers;

	tasksData.reserve(numBins);
	for (uint32_t i = 0; i < numBins; ++i)
		tasksData.push_back(TaskData{ i });

	// sorting is performed in the following manner: first biggest bins are dumped, then smaller; but the sorting is performed basing on the biggest sample only
	std::sort(tasksData.begin(), tasksData.end(), [&](const TaskData& a, const TaskData& b) { return samplesBeginSize[a.binId] > samplesBeginSize[b.binId]; });
}

template<unsigned SIZE>
void Merger<SIZE>::serializeNormalizationAndSave()
{
	typedef StatisticsParams::NormalizationMethod NormalizationMethod;

	std::vector<uint8_t> deseq2NormalizationData, frequencyNormalizationData, quantileNormalizationData;

	if (params.statisticsParams.normalizationMethod == NormalizationMethod::deseq2 || params.statisticsParams.learnDeseq2)
		normalizationLearning.serialize(NormalizationMethod::deseq2, deseq2NormalizationData);
	normalizationLearning.serialize(NormalizationMethod::frequency_count, frequencyNormalizationData);
	normalizationLearning.serialize(NormalizationMethod::quantile, quantileNormalizationData);

	MatrixStatsWriter stats_writer(params.mkmcParams.normLearningBinFile);
	if (params.statisticsParams.normalizationMethod == StatisticsParams::NormalizationMethod::deseq2 || params.statisticsParams.learnDeseq2)
		stats_writer.Add(StatisticsParams::getNormalizationMethodStreamName(NormalizationMethod::deseq2), deseq2NormalizationData);
	stats_writer.Add(StatisticsParams::getNormalizationMethodStreamName(NormalizationMethod::frequency_count), frequencyNormalizationData);
	stats_writer.Add(StatisticsParams::getNormalizationMethodStreamName(NormalizationMethod::quantile), quantileNormalizationData);
}

template<unsigned SIZE>
template<typename PerformGenerate_T>
void Merger<SIZE>::callThreads(const kmcdb::Config& config, std::vector<std::vector<uint64_t>>& tot_cnts, uint64_t totKmersAllSamples)
{
	std::vector<std::thread> threads(params.mkmcParams.nThreads);

	PerformGenerate_T::initWriters(params, config);
	progress_bar = std::make_unique<ProgressBar>(params.mkmcParams.verbosity_level == 0 ? 0 : totKmersAllSamples, "Merging", std::cerr, params.mkmcParams.verbosity_level == 0);

	for (uint32_t i_thred = 0; i_thred < params.mkmcParams.nThreads; ++i_thred)
		threads[i_thred] = std::thread([this, &tot_cnts_thread = tot_cnts[i_thred]]
			{ this->operator() < PerformGenerate_T > (tot_cnts_thread); });

	for (std::thread& thread : threads)
		thread.join();

	PerformGenerate_T::closeWriters();
}

template<unsigned SIZE>
bool Merger<SIZE>::mergeParallel()
{
	uint64_t totKmersAllSamples = 0;
	size_t biggestSample = 0;
	openReadersAndVerifySamples(totKmersAllSamples, biggestSample); // out arguments
	fillTaskData(biggestSample);

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

	std::vector<std::string> sampleNames;
	sampleNames.reserve(params.mkmcParams.samples.size());
	for (const auto& sample : params.mkmcParams.samples)
		sampleNames.push_back(sample.name);

	std::vector<std::vector<uint64_t>> tot_cnts(params.mkmcParams.nThreads,
		std::vector<uint64_t>(sampleNames.size()));

	auto doCallThreads = [&]<typename... Generators_T>() -> void
	{
		using PerformGenerate_T = PerformGenerate<BinFileGenerator, Generators_T...>;
		callThreads<PerformGenerate_T>(config, tot_cnts, totKmersAllSamples);
	};

	GeneratorVecToTemplate<>::callTemplateFunction(params.mkmcParams.outputFileTypes, doCallThreads);

	if (params.mkmcParams.totCntGeneration)
		StoreTotCnt(tot_cnts, sampleNames, params);

	if (!matrixNotEmpty)
	{
		if (params.filterParams.filterKmersSequences)
			Logger::Inst().Log("Error: output matrix is empty; try to relax --ci, --cx, --thr, or --thr_rat parameters, --flt file contains too few k-mers, or some samples are too small.");
		else
			Logger::Inst().Log("Error: output matrix is empty; try to relax --ci, --cx, --thr, or --thr_rat parameters, or some samples are too small.");

		Logger::Inst().Log("Info: removing unnecessary output files.", 2);
		// keep params.mkmcParams.outputFileTotCnt, as its contents are reasonable
		std::filesystem::remove(params.mkmcParams.normLearningBinFile);
		std::filesystem::remove(params.mkmcParams.outputMatrixBinFile);
		std::filesystem::remove(params.mkmcParams.outputMatrixFile);
		std::filesystem::remove(params.mkmcParams.outputFASTAFile);

		return false;
	}

	serializeNormalizationAndSave();
	return true;
}


template<unsigned SIZE>
template<typename PerformGenerate_T>
void Merger<SIZE>::operator()(std::vector<uint64_t>& tot_cnts)
{
	PerformGenerate_T performGenerate(params);
	TaskData taskData;
	using ParameterizedKmersSamplesStruct = KmersSamplesStruct<SIZE>;

	while (tasksPool.getTask(taskData))
	{
		StatisticsParams::NormalizationLearning currentBinNormalizationLearning;
		if (params.statisticsParams.normalizationMethod == StatisticsParams::NormalizationMethod::deseq2 || params.statisticsParams.learnDeseq2)
			currentBinNormalizationLearning.register_method(StatisticsParams::NormalizationMethod::deseq2);
		currentBinNormalizationLearning.register_method(StatisticsParams::NormalizationMethod::frequency_count);
		currentBinNormalizationLearning.register_method(StatisticsParams::NormalizationMethod::quantile);
		currentBinNormalizationLearning.set_no_series(params.mkmcParams.samples.size());
		currentBinNormalizationLearning.initialize();

		performGenerate.setBinId(taskData.binId);

		bool matrixNotEmptyInBin = false;
		if (params.filterParams.filterKmersSequences)
		{
			assert(sequencesToFilterReader);
			using Filters = PerformFilter<FilterCountThreshold<ParameterizedKmersSamplesStruct>, FilterSequences<ParameterizedKmersSamplesStruct>>;
			matrixNotEmptyInBin = mergeToGenerators<PerformGenerate_T, Filters>(taskData.binId, tot_cnts, currentBinNormalizationLearning, performGenerate, sequencesToFilterReader->GetBin(taskData.binId));
		}
		else
		{
			using Filters = PerformFilter<FilterCountThreshold<ParameterizedKmersSamplesStruct>>;
			matrixNotEmptyInBin = mergeToGenerators<PerformGenerate_T, Filters>(taskData.binId, tot_cnts, currentBinNormalizationLearning, performGenerate, nullptr);
		}
		normalizationLearningMutex.lock();
		normalizationLearning.merge_with(&currentBinNormalizationLearning, &currentBinNormalizationLearning + 1);
		matrixNotEmpty |= matrixNotEmptyInBin;
		normalizationLearningMutex.unlock();
	}
}
