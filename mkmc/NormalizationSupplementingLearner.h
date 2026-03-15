#pragma once

#include "parameters.h"
#include "Logger.h"
#include "TextFileWritingUtilities.h"
#include <thread>
#include <vector>
#include <algorithm>



template<unsigned SIZE>
class NormalizationSupplementingLearner
{
	const Params& params;
	const std::vector<StatisticsParams::NormalizationMethod>& methodsToLearn;

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

	std::vector<std::string> samples;

public:
	NormalizationSupplementingLearner(const Params& params, const std::vector<StatisticsParams::NormalizationMethod>& methodsToLearn) :
		params(params),
		methodsToLearn(methodsToLearn),
		tasksPool(tasksData)
	{
		for (auto method : methodsToLearn)
			normalizationLearning.register_method(method);
		normalizationLearning.set_no_series(params.mkmcParams.samples.size());
		normalizationLearning.initialize();
	}

	void learnDeseq2Parallel();

private:
	std::unique_ptr<kmcdb::MetadataReader> matrixMetadataReader;
	std::unique_ptr<kmcdb::ReaderSortedPlainForListing<uint64_t>> matrixReader;

	void openReaders();
	void operator()();
	void fillTaskData();
	void serializeNormalizationAndSave();
};



class NormalizationSupplementingLearnerRunner
{
	Params& params;
	const std::vector<StatisticsParams::NormalizationMethod>& methodsToLearn;
public:
	NormalizationSupplementingLearnerRunner(Params& params, const std::vector<StatisticsParams::NormalizationMethod>& methodsToLearn) :
		params(params),
		methodsToLearn(methodsToLearn)
	{}
	template<unsigned SIZE>
	void Run()
	{
		NormalizationSupplementingLearner<SIZE> deseq2LearnRunner(params, methodsToLearn);
		deseq2LearnRunner.learnDeseq2Parallel();
	}
};



template<unsigned SIZE>
void NormalizationSupplementingLearner<SIZE>::learnDeseq2Parallel()
{
	openReaders();
	fillTaskData();

	matrixReader->GetSampleNames(samples);

	std::vector<std::thread> threads(params.mkmcParams.nThreads);

	for (uint32_t i_thred = 0; i_thred < params.mkmcParams.nThreads; ++i_thred)
		threads[i_thred] = std::thread([this] { this->operator() (); });

	for (std::thread& thread : threads)
		thread.join();

	serializeNormalizationAndSave();
}



template<unsigned SIZE>
void NormalizationSupplementingLearner<SIZE>::openReaders()
{
	try
	{
		matrixMetadataReader = std::make_unique<kmcdb::MetadataReader>(params.mkmcParams.outputMatrixBinFile, false);
		matrixReader = std::make_unique<kmcdb::ReaderSortedPlainForListing<uint64_t>>(*matrixMetadataReader);
	}
	catch (const std::runtime_error& ex)
	{
		Logger::Inst().Log(std::string("Error: ") + ex.what());
		exit(1);
	}
}


template<unsigned SIZE>
void NormalizationSupplementingLearner<SIZE>:: operator()()
{
	TaskData taskData;

	while (tasksPool.getTask(taskData))
	{
		StatisticsParams::NormalizationLearning currentBinNormalizationLearning;
		for (auto method : methodsToLearn)
			currentBinNormalizationLearning.register_method(method);
		currentBinNormalizationLearning.set_no_series(params.mkmcParams.samples.size());
		currentBinNormalizationLearning.initialize();

		// iterate over k-mers in the current bin
		{
			std::vector<uint64_t> kMersCounts(samples.size());
			kMersCounts.resize(samples.size());

			auto bin = matrixReader->GetBin(taskData.binId);

			kmcdb::CKmer<SIZE> kmer;
			while (bin->NextKmer(kmer, kMersCounts.data()))
				currentBinNormalizationLearning.add_entry(kMersCounts);
		}

		normalizationLearningMutex.lock();
		normalizationLearning.merge_with(&currentBinNormalizationLearning, &currentBinNormalizationLearning + 1);
		normalizationLearningMutex.unlock();
	}
}



template<unsigned SIZE>
void NormalizationSupplementingLearner<SIZE>::fillTaskData()
{
	tasksData.reserve(matrixMetadataReader->GetConfig().num_bins);
	for (uint32_t i = 0; i < matrixMetadataReader->GetConfig().num_bins; ++i)
		tasksData.push_back(TaskData{ i });

	// currently, in contrast to Merger, tasks are not sorted basing on problem size
}

template<unsigned SIZE>
void NormalizationSupplementingLearner<SIZE>::serializeNormalizationAndSave()
{
	// Prepared for multiple normalization methods learning. Creates (or replaces) supplementary file.

	typedef StatisticsParams::NormalizationMethod NormalizationMethod;

	Logger::Inst().Log("Info: generating temporary file with data for normalization " + params.mkmcParams.normLearningBinFileSupplemented + ", supplemented with DESeq2 data.", 2);

	auto methodsMaybyPreviouslyLearned = StatisticsParams::getAllSupportedNormalizationMethods();
	auto alwaysLearnedNormalizationMethods = StatisticsParams::getAlwaysLearnedNormalizationMethods();

	for (auto method : alwaysLearnedNormalizationMethods)
		methodsMaybyPreviouslyLearned.erase(std::find(methodsMaybyPreviouslyLearned.begin(), methodsMaybyPreviouslyLearned.end(), method));
	for (auto method : methodsToLearn)
		methodsMaybyPreviouslyLearned.erase(std::find(methodsMaybyPreviouslyLearned.begin(), methodsMaybyPreviouslyLearned.end(), method));

	std::vector<uint8_t> normalizationData;

	MatrixStatsWriter stats_writer(params.mkmcParams.normLearningBinFileSupplemented);

	MatrixStatsReader stats_reader(params.mkmcParams.normLearningBinFile);
	bool success = true;
	for (auto method : alwaysLearnedNormalizationMethods)
	{
		success &= stats_reader.Get(StatisticsParams::getNormalizationMethodStreamName(method), normalizationData);
		stats_writer.Add(StatisticsParams::getNormalizationMethodStreamName(method), normalizationData);
	}

	for (auto method : methodsToLearn)
	{
		assert(std::find(alwaysLearnedNormalizationMethods.begin(), alwaysLearnedNormalizationMethods.end(), method) == alwaysLearnedNormalizationMethods.end());
		// assert(!stats_reader.Get(StatisticsParams::getNormalizationMethodStreamName(method), normalizationData)); - possible if --learn-q or --learn-deseq are used
		normalizationLearning.serialize(method, normalizationData);
		stats_writer.Add(StatisticsParams::getNormalizationMethodStreamName(method), normalizationData);
	}

	if (!success)
	{
		Logger::Inst().Log("Error: cannot read normalization data for supplementing normalization data.");
		std::filesystem::remove(params.mkmcParams.normLearningBinFileSupplemented);
		exit(1);
	}

	// read previously learned stats from input binary file
	for (auto method : methodsMaybyPreviouslyLearned)
	{
		// We do not check success, because acutally we do not know, if this method was utulized
		if (stats_reader.Get(StatisticsParams::getNormalizationMethodStreamName(method), normalizationData))
			stats_writer.Add(StatisticsParams::getNormalizationMethodStreamName(method), normalizationData);
	}
}
