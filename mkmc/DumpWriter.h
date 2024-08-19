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

	void StoreHeader(const std::vector<std::string>& sample_names, const std::string& first_item = "k-mer");

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

//collection of different methods to store k-mer and values
//mkokot_TODO: consider moving to different file?
//It seems to make it work without the need to specify VALUE_T directly when calling OutputBuffer::StoreKmer we need to use generic lambdas
namespace StoreMethods
{
	namespace detail
	{
		template<typename VALUE_T>
		size_t store_single_value(const VALUE_T& val, char*& out, char term)
		{
			size_t r{};
			if constexpr (std::is_integral_v<VALUE_T>)
				r = refresh::int_to_pchar(val, out, term);
			else if constexpr (std::is_floating_point_v<VALUE_T>)
				r = refresh::real_to_pchar(val, out, 6, term);
			else
				static_assert(!sizeof(VALUE_T), "Unsupported type");

			out += r;
			return r;
		}

		inline size_t store_kmer(const std::string& kmerSeq, char* &out, char term)
		{
			std::memcpy(out, kmerSeq.data(), kmerSeq.length());
			out += kmerSeq.length();
			*out++ = term;

			return kmerSeq.length() + 1;
		}
	}
	inline auto AsMatrixRow = []<typename VALUE_T>(const std::string& kmerSeq, const std::vector<VALUE_T>& cnts, char* out) -> size_t
	{
		size_t res = detail::store_kmer(kmerSeq, out, '\t');

		for (size_t i = 0; i < cnts.size() - 1; ++i)
			res += detail::store_single_value(cnts[i], out, '\t');
		res += detail::store_single_value(cnts.back(), out, '\n');

		return res;
	};

	inline auto AsMatrixRow_single_val = []<typename VALUE_T>(const std::string & kmerSeq, const VALUE_T cnt, char* out) -> size_t
	{
		size_t res = detail::store_kmer(kmerSeq, out, '\t');

		res += detail::store_single_value(cnt, out, '\n');

		return res;
	};

	inline auto AsFastaRecord = []<typename VALUE_T>(const std::string& kmerSeq, const std::vector<VALUE_T>&cnts, char* out) -> size_t
	{
		out[0] = '>';
		out[1] = '\n';
		size_t res = 2;
		out += 2;

		res += detail::store_kmer(kmerSeq, out, '\n');

		return res;
	};
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

	template<typename VALUE_T, typename STORE_METHOD>
	void StoreKmer(const std::string& kmerSeq, const std::vector<VALUE_T>& cnts, const STORE_METHOD& storeMethod)
	{
		if (out_buff_pos + max_line_len > buff_owner.size())
		{
			writer.Write(buff, out_buff_pos);
			out_buff_pos = 0;
		}
		out_buff_pos += storeMethod(kmerSeq, cnts, buff + out_buff_pos);
	}

	template<typename VALUE_T, typename STORE_METHOD>
	void StoreKmer(const std::string& kmerSeq, const VALUE_T cnt, const STORE_METHOD& storeMethod)
	{
		if (out_buff_pos + max_line_len > buff_owner.size())
		{
			writer.Write(buff, out_buff_pos);
			out_buff_pos = 0;
		}
		out_buff_pos += storeMethod(kmerSeq, cnt, buff + out_buff_pos);
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
