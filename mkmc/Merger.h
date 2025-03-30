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

	DumpWriter writer(params.mkmcParams.outputFileTotCnt, false);
	writer.StoreHeader(sample_names, "measure");

	std::string first_col = "tot_cnt";
	auto first_col_len = first_col.length();

	auto max_line_len = first_col_len + 1 + params.mkmcParams.samples.size() * (refresh::numeric_conversion_max_length<uint64_t>() + 1);
	OutputBuffer out(writer, max_line_len);
	out.StoreKmer(first_col, tot_cnts[0], StoreMethods::AsMatrixRow);
}

template<unsigned SIZE>
class Merger
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
	void serializeNormalizationAndSave();

	template<typename PerformGenerate_T>
	void callThreads(const kmcdb::Config& config, std::vector<std::vector<uint64_t>>& tot_cnts);

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
	void mergeToGenerators(uint32_t binId, std::vector<uint64_t>& tot_cnts, StatisticsParams::NormalizationLearning& currentBinNormalizationLearning, PerformGenerate_T& fileGenerators, kmcdb::BinReaderSortedWithLUTForListing<uint64_t>* bin);

public:
	Merger(const Params& params) :
		params(params), tasksPool(tasksData)
	{
		if (params.statisticsParams.normalizationMethod == StatisticsParams::NormalizationMethod::deseq2)
			normalizationLearning.register_method(StatisticsParams::NormalizationMethod::deseq2);
		normalizationLearning.register_method(StatisticsParams::NormalizationMethod::frequency_count);
		normalizationLearning.register_method(StatisticsParams::NormalizationMethod::quantile);
		normalizationLearning.set_no_series(params.mkmcParams.samples.size());
		normalizationLearning.initialize();
	}

	void mergeParallel();

	template<typename PerformGenerate_T>
	void operator()(std::vector<uint64_t>& tot_cnts);
};



template<unsigned SIZE>
template<typename PerformGenerate_T, typename Filters_T>
void Merger<SIZE>::mergeToGenerators(uint32_t binId, std::vector<uint64_t>& tot_cnts, StatisticsParams::NormalizationLearning& normalizationLearning, PerformGenerate_T& fileGenerators, kmcdb::BinReaderSortedWithLUTForListing<uint64_t>* bin)
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
						for (size_t i = 0; i < kMersCounts.size(); ++i)
							tot_cnts[i] += kMersCounts[i];

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
		for (size_t i = 0; i < kMersCounts.size(); ++i)
			tot_cnts[i] += kMersCounts[i];

		fileGenerators.writeKmer(KmersSamplesStruct<SIZE>{ minKmer, kMersCounts });
		normalizationLearning.add_entry(kMersCounts);
	}
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
inline void Merger<SIZE>::fillTaskData()
{
	uint64_t biggestSampleKmersCount = 0;

	std::vector<uint64_t> samplesBeginSize;

	//mkokot_TODO: Macku ten reserve to chyba na l. sampli powiniene byc a nie na liczbe binow, prawda?
	//Generalnie na GetNBins nie mozna polegac bo kmc moglo sobie wybrac inna liczbe binow (tzn. 1 jak wlacza small k opt)
	//samplesMetadata.reserve(params.stage1Params.GetNBins());
	//samplesReaders.reserve(params.stage1Params.GetNBins());
	samplesMetadata.reserve(params.mkmcParams.kmcOutputFiles.size());
	samplesReaders.reserve(params.mkmcParams.kmcOutputFiles.size());

	uint64_t totKmersAllSamples = 0;

	for (size_t i = 0; i < params.mkmcParams.kmcOutputFiles.size(); ++i)
	{
		try
		{
			samplesMetadata.emplace_back(std::make_unique<kmcdb::MetadataReader>(params.mkmcParams.kmcOutputFiles[i], true));
			kmcdb::MetadataReader& metadata_reader = *samplesMetadata.back();
			samplesReaders.emplace_back(std::make_unique<kmcdb::ReaderSortedWithLUTForListing<uint64_t>>(metadata_reader));
			kmcdb::ReaderSortedWithLUTForListing<uint64_t>& reader = *samplesReaders.back();

			if (i == 0) //we need to get number of bins from readers, not from config because KMC could enable small k opt and use 1 bin
				samplesBeginSize.resize(metadata_reader.GetConfig().num_bins);

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

	tasksData.reserve(samplesMetadata.front()->GetConfig().num_bins);
	for (uint32_t i = 0; i < samplesMetadata.front()->GetConfig().num_bins; ++i)
	{
		tasksData.push_back(TaskData{ i });
	}

	// sorting is performed in the following manner: first biggest bins are dumped, then smaller; but the sorting is performed basing on the biggest sample only
	std::sort(tasksData.begin(), tasksData.end(), [&](const TaskData& a, const TaskData& b) { return samplesBeginSize[a.binId] > samplesBeginSize[b.binId]; });
}

template<unsigned SIZE>
void Merger<SIZE>::serializeNormalizationAndSave()
{
	std::vector<uint8_t> deseq2NormalizationData, frequencyNormalizationData, quantileNormalizationData;

	if (params.statisticsParams.normalizationMethod == StatisticsParams::NormalizationMethod::deseq2)
		normalizationLearning.serialize(StatisticsParams::NormalizationMethod::deseq2, deseq2NormalizationData);
	normalizationLearning.serialize(StatisticsParams::NormalizationMethod::frequency_count, frequencyNormalizationData);
	normalizationLearning.serialize(StatisticsParams::NormalizationMethod::quantile, quantileNormalizationData);

	MatrixStatsWriter stats_writer(params.mkmcParams.normLearningBinFile);
	if (params.statisticsParams.normalizationMethod == StatisticsParams::NormalizationMethod::deseq2)
		stats_writer.Add(params.statisticsParams.normDeseq2StreamName, deseq2NormalizationData);
	stats_writer.Add(params.statisticsParams.normFrequencyStreamName, frequencyNormalizationData);
	stats_writer.Add(params.statisticsParams.normQuantileStreamName, quantileNormalizationData);
}

template<unsigned SIZE>
template<typename PerformGenerate_T>
void Merger<SIZE>::callThreads(const kmcdb::Config& config, std::vector<std::vector<uint64_t>>& tot_cnts)
{
	std::vector<std::thread> threads(params.mkmcParams.nThreads);

	PerformGenerate_T::initWriters(params, config);
	for (uint32_t i_thred = 0; i_thred < params.mkmcParams.nThreads; ++i_thred)
		threads[i_thred] = std::thread([this, &tot_cnts_thread = tot_cnts[i_thred]]
			{ this->operator() < PerformGenerate_T > (tot_cnts_thread); });

	for (std::thread& thread : threads)
		thread.join();

	PerformGenerate_T::closeWriters();
}

template<unsigned SIZE>
void Merger<SIZE>::mergeParallel()
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

	std::vector<std::string> sampleNames;
	sampleNames.reserve(params.mkmcParams.samples.size());
	for (const auto& sample : params.mkmcParams.samples)
		sampleNames.push_back(sample.name);

	std::vector<std::vector<uint64_t>> tot_cnts(params.mkmcParams.nThreads,
		std::vector<uint64_t>(sampleNames.size()));

	auto doCallThreads = [&]<typename... Generators_T>() -> void
	{
		using PerformGenerate_T = PerformGenerate<BinFileGenerator, Generators_T...>;
		callThreads<PerformGenerate_T>(config, tot_cnts);
	};

	GeneratorVecToTemplate<>::callTemplateFunction(params.mkmcParams.outputFileTypes, doCallThreads);

	if (params.mkmcParams.totCntGeneration)
		StoreTotCnt(tot_cnts, sampleNames, params);

	if (params.statisticsParams.generateNormalization || params.statisticsParams.generateEntropy || !params.statisticsParams.classificationMethods.empty())
		serializeNormalizationAndSave();
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
		if (params.statisticsParams.normalizationMethod == StatisticsParams::NormalizationMethod::deseq2)
			currentBinNormalizationLearning.register_method(StatisticsParams::NormalizationMethod::deseq2);
		currentBinNormalizationLearning.register_method(StatisticsParams::NormalizationMethod::frequency_count);
		currentBinNormalizationLearning.register_method(StatisticsParams::NormalizationMethod::quantile);
		currentBinNormalizationLearning.set_no_series(params.mkmcParams.samples.size());
		currentBinNormalizationLearning.initialize();

		performGenerate.setBinId(taskData.binId);

		if (params.filterParams.filterKmersSequences)
		{
			assert(sequencesToFilterReader);
			using Filters = PerformFilter<FilterCountThreshold<ParameterizedKmersSamplesStruct>, FilterSequences<ParameterizedKmersSamplesStruct>>;
			mergeToGenerators<PerformGenerate_T, Filters>(taskData.binId, tot_cnts, currentBinNormalizationLearning, performGenerate, sequencesToFilterReader->GetBin(taskData.binId));
		}
		else
		{
			using Filters = PerformFilter<FilterCountThreshold<ParameterizedKmersSamplesStruct>>;
			mergeToGenerators<PerformGenerate_T, Filters>(taskData.binId, tot_cnts, currentBinNormalizationLearning, performGenerate, nullptr);
		}
		normalizationLearningMutex.lock();
		normalizationLearning.merge_with(&currentBinNormalizationLearning, &currentBinNormalizationLearning + 1);
		normalizationLearningMutex.unlock();
	}
}
