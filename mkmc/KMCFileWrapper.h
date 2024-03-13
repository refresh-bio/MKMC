#pragma once

#include <string>
#include <cstdint>
#include "../kmc/kmc_api/kmc_file.h"
#include "Filter.h"



class KMCFileWrapper
{
public:
	using kmer_t = CKmerAPI;

private:
	std::unique_ptr<CKMCFile> kmc_file;
	size_t tot_kmers;
	size_t cur_kmer_no{};
	kmer_t cur;
	size_t cur_count;
	uint32_t k;
public:
	KMCFileWrapper(KMCFileWrapper&&) = default;
	KMCFileWrapper& operator=(KMCFileWrapper&&) = default;

	KMCFileWrapper(const std::string& path);
	uint32_t GetK() const
	{
		return k;
	}
	bool Finished() const
	{
		return cur_kmer_no > tot_kmers;
	}
	const kmer_t& First() const
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
			std::cerr << "Error: critical, this should not happen, details: " << __FILE__ << "(" << __LINE__ << ")." << std::endl;
			exit(1);
		}
	}
	++cur_kmer_no;
}