#include "KMCFileWrapper.h"



KMCFileWrapper::KMCFileWrapper(const std::string& path)
{
	kmc_file = std::make_unique<CKMCFile>();
	if (!kmc_file->OpenForListing(path))
	{
		std::cerr << "Error: cannot open kmc database " << path << "\n";
		exit(1);
	}
	if (kmc_file->IsKMC2() == true)
	{
		std::cerr << "Error: kmc database not sorted: " << path << "\n";
		std::cerr << "Hint: KMC database may be sorted with kmc_tools transform <input> sort <output>\n";
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

inline void KMCFileWrapper::Next()
{
	auto read_cnt = [this] {
#ifdef __APPLE__
		uint64 cnt;
		bool res = kmc_file->ReadNextKmer(cur, cnt);
		cur_count = cnt;
		return res;
#else
		return kmc_file->ReadNextKmer(cur, cur_count);
#endif
	};

	if (!read_cnt())
	{
		// for the last one k-mer (cur_kmer_no == tot_kmers), it is legal to call Next, but ReadNextKmer will fail
		if (cur_kmer_no != tot_kmers)
		{
			std::cerr << "Error: critical, this should not happen, details: " << __FILE__ << "(" << __LINE__ << ")\n";
			exit(1);
		}
	}
	++cur_kmer_no;
}