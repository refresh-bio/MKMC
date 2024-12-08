#pragma once

#include "parameters.h"
#include <string>

class Start
{
	const Params& params;

	bool canCreateFile(const std::string& path);
	bool canCreateFileInPath(const std::string& path);

public:
	Start(const Params& params) :
		params(params)
	{}

	bool verifyFiles();
};