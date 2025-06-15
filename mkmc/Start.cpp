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
		std::cerr << "Error: Cannot create output file: " << params.mkmcParams.outputFilesTemplate << ". If " << params.mkmcParams.outputFilesTemplate << " is a path to another directory - does the directory exist?" << std::endl;
		return false;
	}

	if (!std::filesystem::exists(params.mkmcParams.tmpPath))
	{
		if (!std::filesystem::create_directory(params.mkmcParams.tmpPath))
		{
			std::cerr << "Error: the specified directory " << params.mkmcParams.tmpPath << " does not exist and it cannot be created." << std::endl;
			return false;
		}
		else
		{
			//if (params.mkmcParams.verbosity_level > 0) // disabled to proper warningPrinted handling
			//{
				std::cerr << "Warning: the specified directory " << params.mkmcParams.tmpPath << " does not exist. It will be temporarily created." << std::endl;
				if (params.mkmcParams.keepTmpFiles)
					std::cerr << "Warning: as " << params.mkmcParams.tmpPath << " directory was created by MKC, the temporary files will not be kept (--keep parameter will be ignored)." << std::endl;
				warningPrinted = true;
			//}
			params.mutableParams.tmpDirCreated = true;
		}
	}
	else if (!std::filesystem::is_directory(params.mkmcParams.tmpPath))
	{
		std::cerr << "Error: " << params.mkmcParams.tmpPath << " exists, but is not a directory." << std::endl;
		return false;
	}

	if (!canCreateFileInPath(params.stage1Params.GetTmpPath()))
	{
		std::cerr << "Error: Cannot create file in the specified working directory: " << params.stage1Params.GetTmpPath() << "." << std::endl;
		if (params.mutableParams.tmpDirCreated)
			std::filesystem::remove(params.mkmcParams.tmpPath);
		return false;
	}

	return true;
}
