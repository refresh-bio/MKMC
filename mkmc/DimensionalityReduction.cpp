#include "DimensionalityReduction.h"
#include "DumpWriter.h"

#include <iostream>


void DimensionalityReduction::add(const uint64_t feature_idx, const std::vector<out_kmcdb_value_type>& featureEntry)
{
	assert(featureEntry.size() == samplesNames.size());
	if (runUMAP)
	{
		for (size_t sample_id = 0; sample_id < featureEntry.size(); ++sample_id)
			umap->add(sample_id, feature_idx, featureEntry[sample_id]);
	}
}


void DimensionalityReduction::runAndStore()
{
	if (runUMAP)
	{
		try
		{
			umap->run(params.statisticsParams.umap_dimensions);
		}
		catch (const std::length_error&)
		{
			std::cerr << "Error: Cannot run UMAP. Try to tight filtering criteria" << std::endl;
		}

		const auto& umap_res = umap->result();
		assert(umap_res.front().size() == params.statisticsParams.umap_dimensions);

		DumpWriter writer(params.mkmcParams.outputFileUMAP, false);
		writer.StoreHeader(samplesNames, "Dimension");

		std::string first_col_prefix = "UMAP";
		auto first_col_len = first_col_prefix.length() + 10; // I assume 10 is more then enough to store component/dimension number

		auto max_line_len = first_col_len + 1 + params.mkmcParams.samples.size() * (refresh::numeric_conversion_max_length<out_kmcdb_value_type>() + 1);

		OutputBuffer out(writer, max_line_len);

		std::vector<out_kmcdb_value_type> values(samplesNames.size());
		for (size_t row = 0; row < params.statisticsParams.umap_dimensions; ++row)
		{
			std::string first_col = first_col_prefix + std::to_string(row + 1);
			for (size_t sample_id = 0; sample_id < samplesNames.size(); ++sample_id)
				values[sample_id] = umap_res[sample_id][row];
			out.StoreKmer(first_col, values, StoreMethods::AsMatrixRow);
		}
	}
}
