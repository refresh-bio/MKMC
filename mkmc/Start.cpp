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

void Start::verifyFiles()
{
	if (!canCreateFile(params.mkmcParams.outputMatrixFiles.front()))
	{
		std::cerr << "Error: Cannot create output file: " << params.mkmcParams.outputMatrixFiles.front() << "." << std::endl;
		exit(1);
	}

	if (!std::filesystem::exists(params.mkmcParams.tmpPath))
	{
		if (!std::filesystem::create_directory(params.mkmcParams.tmpPath))
		{
			std::cerr << "Error: the specified directory " << params.mkmcParams.tmpPath << " does not exist and it cannot be created." << std::endl;
			exit(1);
		}
		else
		{
			std::cerr << "Warning: the specified directory " << params.mkmcParams.tmpPath << " does not exist. It will be temporary created." << std::endl;
			if (params.mkmcParams.keepTmpFiles)
				std::cerr << "Warning: the temporary files will not be kept (-keep parameter will be ignored)." << std::endl;
			params.mutableParams.tmpDirCreated = true;
		}
	}
	else if (!std::filesystem::is_directory(params.mkmcParams.tmpPath))
	{
		std::cerr << "Error: " << params.mkmcParams.tmpPath << "exists, but is not a directory." << std::endl;
		exit(1);
	}

	if (!canCreateFileInPath(params.stage1Params.GetTmpPath()))
	{
		std::cerr << "Error: Cannot create file in specified working directory: " << params.stage1Params.GetTmpPath() << "." << std::endl;
		if (params.mutableParams.tmpDirCreated)
			std::filesystem::remove(params.mkmcParams.tmpPath);
		exit(1);
	}
}
