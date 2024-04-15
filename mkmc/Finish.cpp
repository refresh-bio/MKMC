#include "Finish.h"
#include <filesystem>



void Finish::finishProcessing()
{
	if (!params.mkmcParams.keepTmpFiles)
	{
		for (const auto& kmcOutputFile : params.mkmcParams.kmcOutputFiles)
		{
			std::filesystem::remove(kmcOutputFile + ".kmc_pre");
			std::filesystem::remove(kmcOutputFile + ".kmc_suf");
		}
		if (params.filterParams.filterKmersSequences)
		{
			if (params.mutableParams.createdFastaFile)
			{
				std::filesystem::remove(params.mutableParams.kmersSequencesToFilterOut);
			}
			std::filesystem::remove(params.filterParams.kmersSequencesToFilterOutDB + ".kmc_pre");
			std::filesystem::remove(params.filterParams.kmersSequencesToFilterOutDB + ".kmc_suf");
		}
		if (params.statisticsParams.generateStatistics)
		{
			std::filesystem::remove(params.statisticsParams.normFrequencyFileTmp);
			std::filesystem::remove(params.statisticsParams.normQuantileFileTmp);
		}
	}
	if (params.mutableParams.tmpDirCreated)
		std::filesystem::remove_all(params.mkmcParams.tmpPath);
}
