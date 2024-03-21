#pragma once

#include <string>
#include <cstdint>
#include "../kmc/kmc_api/kmc_file.h"
#include "Filter.h"


template<unsigned SIZE>
class KMCFileWrapper
{
public:

private:
	std::unique_ptr<CKMCFile> kmc_file;
	size_t tot_kmers;
	size_t cur_kmer_no{};
	CKmer<SIZE> cur;
	size_t cur_count;
	uint32_t k;
public:
	KMCFileWrapper(KMCFileWrapper&&) = default;
	KMCFileWrapper& operator=(KMCFileWrapper&&) = default;
	KMCFileWrapper(const std::string& path, const std::string &mapStatsFileName, uint32_t binId);

	uint32_t GetK() const
	{
		return k;
	}
	bool Finished() const
	{
		return cur_kmer_no > tot_kmers;
	}
	const CKmer<SIZE>& First() const
	{
		return cur;
	}
	size_t FirstCount() const
	{
		return cur_count;
	}
	size_t GetTotKmers() const
	{
		return tot_kmers;
	}
	void Next();
	~KMCFileWrapper() noexcept
	{
		if (kmc_file)
			kmc_file->Close();
	}
};


template<unsigned SIZE>
KMCFileWrapper<SIZE>::KMCFileWrapper(const std::string& path, const std::string& mapStatsFileName, uint32_t binId)
{
	kmc_file = std::make_unique<CKMCFile>(true);
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

	cur.clear();

	if (!Finished())
		Next();
}


template<unsigned SIZE>
inline void KMCFileWrapper<SIZE>::Next()
{
	auto read_cnt = [this] {
#ifdef __APPLE__
		uint64 cnt;
		bool res = kmc_file->ReadNextKmer(cur, cnt);
		cur_count = cnt;
		return res;
#else
		return kmc_file->ReadNextKmerFromBin(cur, cur_count);
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
