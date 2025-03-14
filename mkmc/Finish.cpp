#include "Finish.h"
#include <filesystem>



void Finish::finishProcessing()
{
	if (!params.mkmcParams.keepTmpFiles)
	{
		for (const auto& kmcOutputFile : params.mkmcParams.kmcOutputFiles)
		{
			std::filesystem::remove(kmcOutputFile);
		}
		if (params.filterParams.filterKmersSequences)
		{
			if (params.mutableParams.createdFastaFile)
			{
				std::filesystem::remove(params.mutableParams.kmersSequencesToFilterOut);
			}
			std::filesystem::remove(params.filterParams.kmersSequencesToFilterOutDB);
		}

		// currently MKMC runs in bulk mode, thus user rather won't need binary files
		std::filesystem::remove(params.mkmcParams.outputMatrixBinFile);
		if (params.statisticsParams.generateNormalization || params.statisticsParams.generateEntropy || !params.statisticsParams.classificationMethods.empty())
			std::filesystem::remove(params.mkmcParams.normLearningBinFile);
		if (params.statisticsParams.normalizationLearningWasSupplemented)
			std::filesystem::remove(params.mkmcParams.normLearningBinFileSupplemented); // will be useful after modularization
	}
	if (params.mutableParams.tmpDirCreated)
	{
		std::filesystem::remove_all(params.mkmcParams.tmpPath);
	}
}
