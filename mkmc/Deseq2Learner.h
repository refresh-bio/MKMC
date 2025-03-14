#pragma once

#include "parameters.h"
#include "DumpWriter.h"
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
		std::cerr << "\nStarting learning for DESeq2 normalization\n";

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
		std::cerr << "Error: " << ex.what() << std::endl;
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
	std::vector<uint8_t> deseq2NormalizationData;

	normalizationLearning.serialize(StatisticsParams::NormalizationMethod::deseq2, deseq2NormalizationData);

	// read another stats from input binary file
	MatrixStatsReader stats_reader(params.mkmcParams.normLearningBinFile);
	std::vector<uint8_t> deseq2DummyNormalizationData, freqNormalizationData, quantileNormalizationData;

	assert(!stats_reader.Get(params.statisticsParams.normDeseq2StreamName, deseq2DummyNormalizationData));
	bool success = stats_reader.Get(params.statisticsParams.normFrequencyStreamName, freqNormalizationData);
	success &= stats_reader.Get(params.statisticsParams.normQuantileStreamName, quantileNormalizationData);

	if (!success)
	{
		std::cerr << "Error: cannot read normalization data for DESeq2\n";
		exit(1);
	}

	MatrixStatsWriter stats_writer(params.mkmcParams.normLearningBinFileSupplemented);
	stats_writer.Add(params.statisticsParams.normDeseq2StreamName, deseq2NormalizationData);
	stats_writer.Add(params.statisticsParams.normFrequencyStreamName, freqNormalizationData);
	stats_writer.Add(params.statisticsParams.normQuantileStreamName, quantileNormalizationData);
}
