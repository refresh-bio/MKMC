#include "KMCFileWrapper.h"



KMCFileWrapper::KMCFileWrapper(const std::string& path)
{
	kmc_file = std::make_unique<CKMCFile>();
	if (!kmc_file->OpenForListing(path))
	{
		std::cerr << "Error: cannot open kmc database " << path << "." << std::endl;
		exit(1);
	}
	if (kmc_file->IsKMC2() == true)
	{
		std::cerr << "Error: kmc database not sorted: " << path << "." << std::endl;
		exit(1);
	}
	CKMCFileInfo kmc_file_info;
	kmc_file->Info(kmc_file_info);
	k = kmc_file_info.kmer_length;
	tot_kmers = kmc_file_info.total_kmers;

	cur = kmer_t(k);

	if (!Finished())
		Next();
}
