#pragma once

#include <vector>
#include <cstdint>
#include "parameters.h"



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
	bool Finished()
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
	void Next();
	~KMCFileWrapper() noexcept
	{
		if (kmc_file)
			kmc_file->Close();
	}
};



class Dump
{
	const Params& params;

	bool allKAreSame(const std::vector<KMCFileWrapper>& samples);
	void openDatabases(std::vector<KMCFileWrapper>& samples);

public:
	Dump(const Params& params) :
		params(params)
	{}

	void dumpToStd();
	void dumpToFile();
};
