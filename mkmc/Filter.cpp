#include "Filter.h"

bool Filter::keepKMer(const std::vector<size_t>& kMersCounts)
{
	std::size_t nAboveThreshold = 0;
	for (std::size_t count : kMersCounts)
	{
		if (count >= params.filterParams.minCountThreshold)
			++nAboveThreshold;
	}
	double fracPresent = static_cast<double>(nAboveThreshold) / kMersCounts.size();

	return fracPresent >= params.filterParams.minKmersAboveThresholdRatio;
}