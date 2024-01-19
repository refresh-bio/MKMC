#pragma once

#include "parameters.h"
#include <vector>



class Filter
{

	const Params& params;

public:
	Filter(const Params& params) :
		params(params)
	{}

	bool keepKMer(std::vector<size_t> kMersCounts);
};