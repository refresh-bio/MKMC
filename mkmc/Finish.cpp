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
	}
	if (params.mutableParams.tmpDirCreated)
		std::filesystem::remove_all(params.mkmcParams.tmpPath);
}
