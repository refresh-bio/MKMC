#include "Finish.h"
#include <filesystem>

void Finish::finishProcessing()
{
	if (!params.mkmcParams.keepTmpFiles)
		for (const auto& toolsOutputFile : params.mkmcParams.kmcOutputFiles)
		{
			std::filesystem::remove(toolsOutputFile + ".kmc_pre");
			std::filesystem::remove(toolsOutputFile + ".kmc_suf");
		}
	if (params.mutableParams.tmpDirCreated)
		std::filesystem::remove_all(params.mkmcParams.tmpPath);
}
