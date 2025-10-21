#include "Finish.h"
#include <filesystem>



void Finish::finishProcessing()
{
	if (!params.mkmcParams.keepTmpFiles)
	{
		Logger::Inst().Log("Info: --keep switch not given, removing temporary files.", 2);
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
		std::filesystem::remove(params.mkmcParams.normLearningBinFile); // sometimes not necessary
		std::filesystem::remove(params.mkmcParams.normLearningBinFileSupplemented); // sometimes not necessary; will be useful after modularization
	}
	else
		Logger::Inst().Log("Info: --keep switch given, keeping temporary files.", 2);
	if (params.mutableParams.tmpDirCreated)
	{
		std::filesystem::remove_all(params.mkmcParams.tmpPath);
		Logger::Inst().Log("Info: removing temporary directory " + params.mkmcParams.tmpPath + " created by MKMC.", 2);
	}
}
