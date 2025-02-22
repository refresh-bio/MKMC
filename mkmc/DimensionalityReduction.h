#pragma once

#include "parameters.h"
#include "lib/refresh/statistics/lib/statistics_umap.h"

#include <vector>
#include <memory>



class DimensionalityReduction {
	Params& params;
	const std::vector<std::string>& samplesNames;

	using out_kmcdb_value_type = double;

	bool runUMAP;

	std::unique_ptr<refresh::umap_direct<out_kmcdb_value_type>> umap;

public:
	DimensionalityReduction(Params& params, const std::vector<std::string>& samplesNames, uint64_t nKmers) :
		params(params),
		samplesNames(samplesNames),
		runUMAP(params.statisticsParams.runUMAP)
	{
		if (runUMAP)
		{
			umap = std::make_unique<refresh::umap_direct<out_kmcdb_value_type>>(samplesNames.size(), nKmers);
			umap->set_params(params.statisticsParams.umap_params);
		}
	}

	void add(const uint64_t feature_idx, const std::vector<out_kmcdb_value_type>& featureEntry);
	void runAndStore();
};
