#include "Filter.h"

bool Filter::keepKMer(const std::vector<size_t>& kMersCounts)
{
	std::size_t nPresent = 0;
	for (std::size_t count : kMersCounts)
	{
		if (count > 0)
			++nPresent;
	}
	double fracPresent = static_cast<double>(nPresent) / kMersCounts.size();

	return fracPresent >= params.mkmcParams.minKmersPresenceThreshold;
}