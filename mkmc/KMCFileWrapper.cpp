#include "KMCFileWrapper.h"



KMCFileWrapper::KMCFileWrapper(const std::string& path, const std::string& mapStatsFileName, uint32_t binId)
{
	kmc_file = std::make_unique<CKMCFile>();
	if (!kmc_file->OpenForListingWithBinOrder(path, mapStatsFileName))
	{
		std::cerr << "Error: cannot open kmc database " << path << "\n";
		exit(1);
	}
	kmc_file->StartBin(binId);
	if (!kmc_file->IsKMC2())
	{
		std::cerr << "Error: this version requires KMC 2 database format: " << path << "\n";
		exit(1);
	}
	CKMCFileInfo kmc_file_info;
	kmc_file->Info(kmc_file_info);
	k = kmc_file_info.kmer_length;

	kmc_file->GetNKmers(binId, tot_kmers);

	cur = kmer_t(k);

	if (!Finished())
		Next();
}

