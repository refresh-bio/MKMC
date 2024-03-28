#pragma once

#include "parameters.h"
#include "KMCFileWrapper.h"
#include "../kmc/kmc_api/kmc_file.h"



template <typename KmersSamplesData>
class FilterKeepAll
{
public:
	FilterKeepAll(const Params& params, uint32_t binId) {}
	bool keepKMer(const KmersSamplesData& kmersData) {
		return true;
	}
};



template <typename KmersSamplesData, typename NextFilter = FilterKeepAll<KmersSamplesData>> 
class FilterCountThreshold
{
	const Params& params;
	NextFilter nextFilter;

public:
	FilterCountThreshold(const Params& params, uint32_t binId) :
		params(params),
		nextFilter(params, binId)
	{}

	bool keepKMer(const KmersSamplesData& kmersData);
};



template <typename KmersSamplesData, typename NextFilter = FilterKeepAll<KmersSamplesData>>
class FilterSequences
{
	const Params& params;
	NextFilter nextFilter;
	KMCFileWrapper<KmersSamplesData::SIZE> kmcFile;

public:
	FilterSequences(const Params& params, uint32_t binId) :
		params(params),
		nextFilter(params, binId),
		kmcFile(params.filterParams.kmersSequencesToFilterOutDB, binId)
	{}

	bool keepKMer(const KmersSamplesData& kmersData);
};



template <typename KmersSamplesData, typename NextFilter>
bool FilterCountThreshold<KmersSamplesData, NextFilter>::keepKMer(const KmersSamplesData& kmersData)
{
	std::size_t nAboveThreshold = 0;
	for (std::size_t count : kmersData.kMersCounts)
	{
		if (count >= params.filterParams.minCountThreshold)
			++nAboveThreshold;
	}
	double fracPresent = static_cast<double>(nAboveThreshold) / kmersData.kMersCounts.size();

	if (fracPresent >= params.filterParams.minKmersAboveThresholdRatio) {
		return nextFilter.keepKMer(kmersData);
	}
	return false;
}



template<typename KmersSamplesData, typename NextFilter>
bool FilterSequences<KmersSamplesData, NextFilter>::keepKMer(const KmersSamplesData& kmersData)
{
	if (kmcFile.Finished())
	{
		return nextFilter.keepKMer(kmersData);
	}

	//assert(!(kmcFile.First() < kmersData.minKmer));
	if (kmersData.minKmer < kmcFile.First())
	{
		return nextFilter.keepKMer(kmersData);
	}
	kmcFile.Next();
	return false;
}
