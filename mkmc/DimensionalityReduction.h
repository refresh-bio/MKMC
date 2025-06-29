#pragma once

#include "parameters.h"
#include "refresh/statistics/lib/statistics_umap.h"
#include "refresh/statistics/lib/statistics_pca.h"

#include <vector>
#include <memory>



class DimensionalityReduction {
	Params& params;
	const std::vector<std::string>& samplesNames;

	using out_kmcdb_value_type = double;

	bool runUMAP;
	bool runPCA;

	std::unique_ptr<refresh::umap_direct<out_kmcdb_value_type>> umap;

	std::unique_ptr<refresh::pca_parallel_add<out_kmcdb_value_type>> pca;

	void store(const std::string& fileName, const std::string& firstColPrefix, const std::vector<std::vector< out_kmcdb_value_type>>& results);
	void store(const std::string& fileName, const std::vector<std::string>& firstCol, const std::vector<std::string>& header, const std::vector<std::vector<out_kmcdb_value_type>>& results);

	void runAndStoreUMAP();
	void runAndStorePCA();

public:
	DimensionalityReduction(Params& params, const std::vector<std::string>& samplesNames, uint64_t nKmers) :
		params(params),
		samplesNames(samplesNames),
		runUMAP(params.statisticsParams.runUMAP),
		runPCA(params.statisticsParams.runPCA)
	{
		if (runUMAP)
		{
			umap = std::make_unique<refresh::umap_direct<out_kmcdb_value_type>>(samplesNames.size(), nKmers);
			umap->set_params(params.statisticsParams.umap_params);
		}
		if (runPCA)
		{
			pca = std::make_unique<refresh::pca_parallel_add<out_kmcdb_value_type>>(samplesNames.size(), nKmers);
		}
	}

	void add(const uint64_t feature_idx, const std::vector<out_kmcdb_value_type>& featureEntry);
	void runAndStore();
};
