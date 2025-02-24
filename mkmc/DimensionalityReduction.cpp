#include "DimensionalityReduction.h"
#include "DumpWriter.h"

#include <iostream>


void DimensionalityReduction::store(const std::string& fileName, const std::string& firstColPrefix, const std::vector<std::vector<out_kmcdb_value_type>>& results)
{
	DumpWriter writer(fileName, false);
	writer.StoreHeader(samplesNames, "Dimension");

	auto first_col_len = firstColPrefix.length() + 10; // I assume 10 is more then enough to store component/dimension number

	auto max_line_len = first_col_len + 1 + params.mkmcParams.samples.size() * (refresh::numeric_conversion_max_length<out_kmcdb_value_type>() + 1);

	OutputBuffer out(writer, max_line_len);

	std::vector<out_kmcdb_value_type> values(samplesNames.size());
	for (size_t row = 0; row < params.statisticsParams.nDimensionReduction; ++row)
	{
		std::string first_col = firstColPrefix + std::to_string(row + 1);
		for (size_t sample_id = 0; sample_id < samplesNames.size(); ++sample_id)
			values[sample_id] = results[sample_id][row];
		out.StoreKmer(first_col, values, StoreMethods::AsMatrixRow);
	}
}


void DimensionalityReduction::runAndStoreUMAP()
{
	try
	{
		umap->run(params.statisticsParams.nDimensionReduction);
	}
	catch (const std::length_error&)
	{
		std::cerr << "Error: Cannot run UMAP. Try to tight filtering criteria" << std::endl;
	}

	const auto& umap_res = umap->result();
	assert(umap_res.front().size() == params.statisticsParams.nDimensionReduction);

	store(params.mkmcParams.outputFileUMAP, "UMAP", umap_res);
}


void DimensionalityReduction::runAndStorePCA()
{
	try
	{
		pca->run(params.statisticsParams.nDimensionReduction);
	}
	catch (...)
	{
		std::cerr << "Error: Cannot run PCA. Try to tight filtering criteria" << std::endl;
	}

	const auto& pca_res = pca->result();
	assert(pca_res.front().size() == params.statisticsParams.nDimensionReduction);

	store(params.mkmcParams.outputFilePCA, "PCA", pca_res);
}


void DimensionalityReduction::add(const uint64_t feature_idx, const std::vector<out_kmcdb_value_type>& featureEntry)
{
	assert(featureEntry.size() == samplesNames.size());
	if (runUMAP)
	{
		for (size_t sample_id = 0; sample_id < featureEntry.size(); ++sample_id)
			umap->add(sample_id, feature_idx, featureEntry[sample_id]);
	}
	if (runPCA)
	{
		PCAAddMutex.lock();
		pca->add_feature(featureEntry.begin(), featureEntry.end());
		PCAAddMutex.unlock();
	}
}


void DimensionalityReduction::runAndStore()
{
	if (params.statisticsParams.runUMAP)
		runAndStoreUMAP();
	if (params.statisticsParams.runPCA)
		runAndStorePCA();
}
