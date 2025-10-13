#include "DimensionalityReduction.h"
#include "TextFileWritingUtilities.h"
#include "Logger.h"

#include <iostream>


void DimensionalityReduction::store(const std::string& fileName, const std::string& firstColPrefix, const std::vector<std::vector<out_kmcdb_value_type>>& results)
{
	TextFileWriter writer(fileName, false);
	writer.StoreHeader(samplesNames, "Dimension");

	MatrixOutputBuffer<out_kmcdb_value_type> out(writer, firstColPrefix.length() + 10, params.mkmcParams.samples.size()); // I assume 10 is more then enough to store component/dimension number

	std::vector<out_kmcdb_value_type> values(samplesNames.size());
	for (size_t row = 0; row < params.statisticsParams.nDimensionReduction; ++row)
	{
		std::string first_col = firstColPrefix + std::to_string(row + 1);
		for (size_t sample_id = 0; sample_id < samplesNames.size(); ++sample_id)
			values[sample_id] = results[sample_id][row];
		out.StoreKmer(first_col, values);
	}
}

void DimensionalityReduction::store(const std::string& fileName, const std::vector<std::string>& firstCol, const std::vector<std::string>& header, const std::vector<std::vector<out_kmcdb_value_type>>& results)
{
	TextFileWriter writer(fileName, false);
	writer.StoreHeader(header, "Statistic");

	size_t firstColLen = 0;
	for (const auto& s : firstCol)
		if (s.length() > firstColLen)
			firstColLen = s.length();

	MatrixOutputBuffer<out_kmcdb_value_type> out(writer, firstColLen, params.mkmcParams.samples.size());

	std::vector<out_kmcdb_value_type> values(samplesNames.size());
	for (size_t row = 0; row < firstCol.size(); ++row)
	{
		for (size_t sample_id = 0; sample_id < samplesNames.size(); ++sample_id)
			values[sample_id] = results[sample_id][row];
		out.StoreKmer(firstCol[row], values);
	}
}

void DimensionalityReduction::runAndStoreUMAP()
{
	bool success = true;
	try
	{
		umap->run(params.statisticsParams.nDimensionReduction);
	}
	catch (const std::length_error&)
	{
		Logger::Inst().Log("Error: Cannot run UMAP. Try to tight filtering criteria.");
		success = false;
	}

	if (success)
	{
		const auto& umap_res = umap->result();
		assert(umap_res.front().size() == params.statisticsParams.nDimensionReduction);

		store(params.mkmcParams.outputFileUMAP, "UMAP", umap_res);
	}
	else
		Logger::Inst().Log("Info: Despite the UMAP failure, MKMC will continue, but no UMAP results will be created.");
}


void DimensionalityReduction::runAndStorePCA()
{
	bool success = true;
	try
	{
		pca->run(params.statisticsParams.nDimensionReduction, params.statisticsParams.pca_mod);
	}
	catch (const std::length_error&)
	{
		Logger::Inst().Log("Error: Cannot run PCA. Try to tight filtering criteria.");
		success = false;
	}
	catch (...)
	{
		Logger::Inst().Log("Error: Unexpected error with PCA running.");
		success = false;
	}

	if (success)
	{
		const auto& pca_res = pca->result();
		assert(pca_res.front().size() == params.statisticsParams.nDimensionReduction);

		store(params.mkmcParams.outputFilePCA, "PCA", pca_res);

		assert(pca->get_explained_variance().size() == pca->get_explained_variance_ratio().size());
		std::vector<std::vector<out_kmcdb_value_type>> variances;
		std::vector<std::string> header;
		for (size_t i = 0; i < samplesNames.size(); ++i)
		{
			header.push_back("D" + std::to_string(i));
			variances.push_back({ pca->get_explained_variance()[i], pca->get_explained_variance_ratio()[i] });
		}

		store(params.mkmcParams.outputFilePCAVariance, std::vector<std::string>{ "variance", "variance_ratio"}, header, variances);
	}
	else
		Logger::Inst().Log("Info: Despite the PCA failure, MKMC will continue, but no PCA results will be created.");
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
		for (size_t sample_id = 0; sample_id < featureEntry.size(); ++sample_id)
			pca->add(sample_id, feature_idx, featureEntry[sample_id]);
	}
}


void DimensionalityReduction::runAndStore()
{
	if (params.statisticsParams.runUMAP)
		runAndStoreUMAP();
	if (params.statisticsParams.runPCA)
		runAndStorePCA();
}
