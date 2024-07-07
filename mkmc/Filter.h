#pragma once

#include "parameters.h"
#include "KMCFileWrapper.h"
#include "../kmc/kmc_api/kmc_file.h"



template <typename KmersSamplesData_T> 
class FilterCountThreshold
{
	const Params& params;

public:
	FilterCountThreshold(const Params& params, kmcdb::BinReaderSortedWithLUTForListing<uint64_t>* /*bin*/) :
		params(params)
	{}

	bool keepKMer(const KmersSamplesData_T& kmersData);
};



template <typename KmersSamplesData_T>
class FilterSequences
{
	const Params& params;
	KMCFileWrapper<KmersSamplesData_T::SIZE> kmcFile;

public:
	FilterSequences(const Params& params, kmcdb::BinReaderSortedWithLUTForListing<uint64_t>* bin) :
		params(params),
		kmcFile(bin)
	{}

	bool keepKMer(const KmersSamplesData_T& kmersData);
};



template<typename Filter_T, typename... NextFilters_T>
class PerformFilter
{
	Filter_T filter;
	PerformFilter<NextFilters_T...> nextPerformFilter;
public:
	PerformFilter(const Params& params, kmcdb::BinReaderSortedWithLUTForListing<uint64_t>* bin) :
		filter(params, bin), nextPerformFilter(params, bin)
	{}

	template<typename KmersSamplesData_T>
	bool keepKMer(const KmersSamplesData_T& kmersData)
	{
		if (!filter.keepKMer(kmersData))
			return false;
		return nextPerformFilter.keepKMer(kmersData);
	}
};



template<typename Filter_T>
class PerformFilter<Filter_T>
{
	Filter_T filter;
public:
	PerformFilter(const Params& params, kmcdb::BinReaderSortedWithLUTForListing<uint64_t>* bin) :
		filter(params, bin)
	{}

	template<typename KmersSamplesData_T>
	bool keepKMer(const KmersSamplesData_T& kmersData)
	{
		return filter.keepKMer(kmersData);
	}
};



template <typename KmersSamplesData_T>
bool FilterCountThreshold<KmersSamplesData_T>::keepKMer(const KmersSamplesData_T& kmersData)
{
	std::size_t nAboveThreshold = 0;
	for (uint64_t count : kmersData.kMersCounts)
	{
		if (count >= params.filterParams.minCountThreshold)
			++nAboveThreshold;
	}

	const double fracPresent = static_cast<double>(nAboveThreshold) / kmersData.kMersCounts.size();
	if (fracPresent >= params.filterParams.minKmersAboveThresholdRatio) {
		return true;
	}
	return false;
}



template<typename KmersSamplesData_T>
bool FilterSequences<KmersSamplesData_T>::keepKMer(const KmersSamplesData_T& kmersData)
{
	if (kmcFile.Finished())
	{
		return false;
	}

	//assert(!(kmcFile.First() < kmersData.minKmer));
	if (kmersData.minKmer < kmcFile.First())
	{
		return false;
	}
	kmcFile.Next();
	return true;
}
