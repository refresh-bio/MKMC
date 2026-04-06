#include "Start.h"
#include <filesystem>
#include <cstdio>



bool Start::canCreateFile(const std::string& path)
{
	FILE* f = fopen(path.c_str(), "wb");
	if (!f)
		return false;
	fclose(f);
	remove(path.c_str());
	return true;
}

bool Start::canCreateFileInPath(const std::string& path)
{
	static const std::string name = "kmc_test.bin"; //Some random name
	if (path.back() == '\\' || path.back() == '/')
		return canCreateFile(path + name);
	else
		return canCreateFile(path + static_cast<char>(std::filesystem::path::preferred_separator) + name);
}

bool Start::verifyFiles(bool& warningPrinted)
{
	if (!canCreateFile(params.mkmcParams.outputFilesTemplate))
	{
		Logger::Inst().Log("Error: cannot create output file: " + params.mkmcParams.outputFilesTemplate + ". If " + params.mkmcParams.outputFilesTemplate + " is a path to another directory - does the directory exist?");
		return false;
	}

	if (!std::filesystem::exists(params.mkmcParams.tmpPath))
	{
		if (!std::filesystem::create_directory(params.mkmcParams.tmpPath))
		{
			Logger::Inst().Log("Error: the specified directory " + params.mkmcParams.tmpPath + " does not exist and it cannot be created.");
			return false;
		}
		else
		{
			Logger::Inst().Log("Warning: the specified directory " + params.mkmcParams.tmpPath + " does not exist. It will be temporarily created.", 1);
			if (params.mkmcParams.keepKMCdbs)
				Logger::Inst().Log("Warning: as " + params.mkmcParams.tmpPath + " directory was created by MKMC, KMC per-sample databases will not be kept (--keep-kmc-temporary-databases parameter).", 1);
			warningPrinted = true;
			params.mutableParams.tmpDirCreated = true;
		}
	}
	else if (!std::filesystem::is_directory(params.mkmcParams.tmpPath))
	{
		Logger::Inst().Log("Error: " + params.mkmcParams.tmpPath + " exists, but is not a directory.");
		return false;
	}

	if (!canCreateFileInPath(params.stage1Params.GetTmpPath()))
	{
		Logger::Inst().Log("Error: cannot create file in the specified working directory: " + params.stage1Params.GetTmpPath() + ".");
		if (params.mutableParams.tmpDirCreated)
			std::filesystem::remove(params.mkmcParams.tmpPath);
		return false;
	}

	return true;
}
