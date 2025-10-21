#pragma once

#include "parameters.h"
#include "Logger.h"
#include "TextFileWritingUtilities.h"
#include <thread>



template<unsigned SIZE>
class Deseq2Learner
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

	std::vector<std::string> samples;

public:
	Deseq2Learner(const Params& params) :
		params(params),
		tasksPool(tasksData)
	{
		normalizationLearning.register_method(StatisticsParams::NormalizationMethod::deseq2);
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



class Deseq2LearnerRunner
{
	Params& params;
public:
	Deseq2LearnerRunner(Params& params) :
		params(params)
	{}
	template<unsigned SIZE>
	void Run()
	{
		Deseq2Learner<SIZE> deseq2LearnRunner(params);
		deseq2LearnRunner.learnDeseq2Parallel();
	}
};



template<unsigned SIZE>
void Deseq2Learner<SIZE>::learnDeseq2Parallel()
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
void Deseq2Learner<SIZE>::openReaders()
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
void Deseq2Learner<SIZE>:: operator()()
{
	TaskData taskData;

	while (tasksPool.getTask(taskData))
	{
		StatisticsParams::NormalizationLearning currentBinNormalizationLearning;
		currentBinNormalizationLearning.register_method(StatisticsParams::NormalizationMethod::deseq2);
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
void Deseq2Learner<SIZE>::fillTaskData()
{
	tasksData.reserve(matrixMetadataReader->GetConfig().num_bins);
	for (uint32_t i = 0; i < matrixMetadataReader->GetConfig().num_bins; ++i)
		tasksData.push_back(TaskData{ i });

	// currently, in contrast to Merger, tasks are not sorted basing on problem size
}

template<unsigned SIZE>
void Deseq2Learner<SIZE>::serializeNormalizationAndSave()
{
	typedef StatisticsParams::NormalizationMethod NormalizationMethod;

	std::vector<uint8_t> deseq2NormalizationData;

	normalizationLearning.serialize(NormalizationMethod::deseq2, deseq2NormalizationData);

	// read another stats from input binary file
	MatrixStatsReader stats_reader(params.mkmcParams.normLearningBinFile);
	std::vector<uint8_t> deseq2DummyNormalizationData, freqNormalizationData, quantileNormalizationData;

	assert(!stats_reader.Get(StatisticsParams::getNormalizationMethodStreamName(NormalizationMethod::deseq2), deseq2DummyNormalizationData));
	bool success = stats_reader.Get(StatisticsParams::getNormalizationMethodStreamName(NormalizationMethod::frequency_count), freqNormalizationData);
	success &= stats_reader.Get(StatisticsParams::getNormalizationMethodStreamName(NormalizationMethod::quantile), quantileNormalizationData);

	if (!success)
	{
		Logger::Inst().Log("Error: cannot read normalization data for DESeq2.");
		exit(1);
	}

	Logger::Inst().Log("Info: generating temporary file with data for normalization " + params.mkmcParams.normLearningBinFileSupplemented + ", supplemented with DESeq2 data.", 2);
	MatrixStatsWriter stats_writer(params.mkmcParams.normLearningBinFileSupplemented);
	stats_writer.Add(StatisticsParams::getNormalizationMethodStreamName(NormalizationMethod::deseq2), deseq2NormalizationData);
	stats_writer.Add(StatisticsParams::getNormalizationMethodStreamName(NormalizationMethod::frequency_count), freqNormalizationData);
	stats_writer.Add(StatisticsParams::getNormalizationMethodStreamName(NormalizationMethod::quantile), quantileNormalizationData);
}
