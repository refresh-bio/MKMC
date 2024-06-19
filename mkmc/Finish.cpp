#include "Finish.h"
#include <filesystem>



void Finish::finishProcessing()
{
	if (!params.mkmcParams.keepTmpFiles)
	{
		for (const auto& kmcOutputFile : params.mkmcParams.kmcOutputFiles)
		{
			std::filesystem::remove(kmcOutputFile);
			std::filesystem::remove(kmcOutputFile);
		}
		if (params.filterParams.filterKmersSequences)
		{
			if (params.mutableParams.createdFastaFile)
			{
				std::filesystem::remove(params.mutableParams.kmersSequencesToFilterOut);
			}
			std::filesystem::remove(params.filterParams.kmersSequencesToFilterOutDB);
			std::filesystem::remove(params.filterParams.kmersSequencesToFilterOutDB);
		}
		if (params.statisticsParams.generateNormalization)
		{
			std::filesystem::remove(params.statisticsParams.normFrequencyFileTmp);
			std::filesystem::remove(params.statisticsParams.normQuantileFileTmp);
			std::filesystem::remove(params.statisticsParams.statsNOutputKmers);
		}
	}
	if (params.mutableParams.tmpDirCreated)
	{
		std::filesystem::remove_all(params.mkmcParams.tmpPath);
	}
}
