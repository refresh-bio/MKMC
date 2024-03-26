#pragma once

#include "parameters.h"
#include <vector>



class FilterKeepAll
{
public:
	FilterKeepAll(const Params& params) {}
	bool keepKMer(const std::vector<size_t>& kMersCounts) {
		return true;
	}
};



template <typename NextFilter = FilterKeepAll>
class FilterCountThreshold
{
	const Params& params;
	NextFilter nextFilter;

public:
	FilterCountThreshold(const Params& params) :
		params(params),
		nextFilter(params)
	{}

	bool keepKMer(const std::vector<size_t>& kMersCounts);
};



template <typename NextFilter>
bool FilterCountThreshold<NextFilter>::keepKMer(const std::vector<size_t>& kMersCounts)
{
	std::size_t nAboveThreshold = 0;
	for (std::size_t count : kMersCounts)
	{
		if (count >= params.filterParams.minCountThreshold)
			++nAboveThreshold;
	}
	double fracPresent = static_cast<double>(nAboveThreshold) / kMersCounts.size();

	if (fracPresent >= params.filterParams.minKmersAboveThresholdRatio) {
		return nextFilter.keepKMer(kMersCounts);
	}
	return false;
}
