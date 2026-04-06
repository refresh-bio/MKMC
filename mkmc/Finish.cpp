#include "Finish.h"
#include <filesystem>



void Finish::finishProcessing()
{
	if (!params.mutableParams.kmcDbsCreated)
	{
		// Do nothing.
	}
	else if (!params.mkmcParams.keepKMCdbs)
	{
		Logger::Inst().Log("Info: --keep-kmc-temporary-databases switch not given, removing KMC databases.", 2);
		for (const auto& kmcOutputFile : params.mkmcParams.kmcOutputFiles)
			std::filesystem::remove(kmcOutputFile); // sometimes not necessary
		if (params.filterParams.filterKmersSequences)
		{
			if (params.mutableParams.createdFastaFile) // We treat a temporary FASTA as KMC DB
				std::filesystem::remove(params.mutableParams.kmersSequencesToFilterOut);
			std::filesystem::remove(params.filterParams.kmersSequencesToFilterOutDB);
		}
	}
	else
		Logger::Inst().Log("Info: --keep-kmc-temporary-databases switch given, keeping KMC databases.", 2);

	if (!params.mkmcParams.reuseDBFiles)
	{
		Logger::Inst().Log("Info: --reuse-db switch not given, removing binary database files.", 2);
		std::filesystem::remove(params.mkmcParams.outputMatrixBinFile);
		std::filesystem::remove(params.mkmcParams.normLearningBinFile); // sometimes not necessary
		std::filesystem::remove(params.mkmcParams.normLearningBinFileSupplemented); // sometimes not necessary
	}
	else
		Logger::Inst().Log("Info: --reuse-db switch given, keeping binary database files.", 2);

	if (params.mutableParams.tmpDirCreated)
	{
		std::filesystem::remove_all(params.mkmcParams.tmpPath);
		Logger::Inst().Log("Info: removing temporary directory " + params.mkmcParams.tmpPath + " created by MKMC.", 2);
	}
}
