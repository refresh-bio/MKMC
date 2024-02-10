#include "Finish.h"
#include <filesystem>

void Finish::finishProcessing()
{
	if (!params.mkmcParams.keepTmpFiles)
		for (const auto& toolsOutputFile : params.mkmcParams.toolsOutputFiles)
		{
			std::filesystem::remove(toolsOutputFile + ".kmc_pre");
			std::filesystem::remove(toolsOutputFile + ".kmc_suf");
		}
}
