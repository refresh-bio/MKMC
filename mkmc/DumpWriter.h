#pragma once

#include <vector>
#include <cstdio>
#include <mutex>
#include <string>
#include <iostream>
#include <algorithm>
#include "kmcdb/kmcdb.h"



class DumpWriter
{
	FILE* out;
	bool use_mutex = false;
	std::mutex mtx; //mkokot_TODO: probably with the queue it would work better, but require more coding

public:
	DumpWriter(const std::string& path, bool use_mutex) :
		out(fopen(path.c_str(), "wb")),
		use_mutex(use_mutex)
	{
		if (!out)
		{
			std::cerr << "Error: cannot open file " << path << "\n";
			exit(1);
		}
		setvbuf(out, nullptr, _IONBF, 0);
	}

	void StoreHeader(const std::vector<std::string>& sample_names);

	void Write(char* data, size_t size)
	{
		if (!use_mutex) {
			fwrite(data, 1, size, out);
			return;
		}
		std::lock_guard<std::mutex> lck(mtx);
		fwrite(data, 1, size, out);
	}

	~DumpWriter()
	{
		if (fclose(out) != 0)
		{
			std::cerr << "Error: some error occurred when closing the output file\n";
			exit(1);
		}
	}
};



class OutputBuffer
{
	DumpWriter& writer;
	size_t max_line_len;
	size_t out_buff_pos{};
	std::vector<char> buff_owner;
	char* buff;

public:
	OutputBuffer(DumpWriter& writer, size_t max_line_len, size_t buff_size = 1ull << 23) :
		writer(writer),
		max_line_len(max_line_len),
		buff_owner((std::max)(buff_size, max_line_len)),
		buff(buff_owner.data())
	{}

	template<unsigned SIZE, typename VALUE_T, typename STORE_METHOD>
	void StoreKmer(const kmcdb::CKmer<SIZE>& kmer, uint64_t kmer_len, const std::vector<VALUE_T>& cnts, const STORE_METHOD& storeMethod)
	{
		if (out_buff_pos + max_line_len > buff_owner.size())
		{
			writer.Write(buff, out_buff_pos);
			out_buff_pos = 0;
		}
		out_buff_pos += storeMethod(kmer, kmer_len, cnts, buff + out_buff_pos);
	}

	template<unsigned SIZE, typename VALUE_T, typename STORE_METHOD>
	void StoreKmer(const kmcdb::CKmer<SIZE>& kmer, uint64_t kmer_len, const VALUE_T cnt, const STORE_METHOD& storeMethod)
	{
		if (out_buff_pos + max_line_len > buff_owner.size())
		{
			writer.Write(buff, out_buff_pos);
			out_buff_pos = 0;
		}
		out_buff_pos += storeMethod(kmer, kmer_len, cnt, buff + out_buff_pos);
	}

	~OutputBuffer()
	{
		if (out_buff_pos)
		{
			writer.Write(buff, out_buff_pos);
			out_buff_pos = 0;
		}
	}
};
