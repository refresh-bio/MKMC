#pragma once

#include <vector>
#include <cstdio>
#include <mutex>
#include <string>
#include <iostream>
#include <algorithm>
#include "kmcdb/kmcdb.h"



class TextFileWriter
{
	FILE* out;
	bool use_mutex = false;
	std::mutex mtx; //mkokot_TODO: probably with the queue it would work better, but require more coding

public:
	TextFileWriter(const std::string& path, bool use_mutex) :
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

	void StoreHeader(const std::vector<std::string>& data_columns_headers, const std::string& first_column_header = "k-mer");

	void Write(char* data, size_t size)
	{
		if (!use_mutex) {
			fwrite(data, 1, size, out);
			return;
		}
		std::lock_guard<std::mutex> lck(mtx);
		fwrite(data, 1, size, out);
	}

	~TextFileWriter()
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
	TextFileWriter& writer;
	size_t max_line_len;
	size_t out_buff_pos{};
	std::vector<char> buff_owner;
	char* buff;

protected:
	char* get_buffer() const
	{
		return buff + out_buff_pos;
	}

	void shift_buffer_after_saving(size_t record_len);

	template<typename VALUE_T>
	size_t store_single_value(const VALUE_T& val, char*& out, char term) const;

	size_t store_kmer(const std::string& kmerSeq, char*& out, char term) const;

	OutputBuffer(TextFileWriter& writer, size_t max_line_len, size_t buff_size = 1ull << 23) :
		writer(writer),
		max_line_len(max_line_len),
		buff_owner((std::max)(buff_size, max_line_len)),
		buff(buff_owner.data())
	{}

	~OutputBuffer();
};

template <typename VALUE_T>
class MatrixOutputBuffer : public OutputBuffer
{
private:
	size_t get_max_record_len(size_t first_col_len, size_t num_columns) const
	{
		//     k-mer                          term(\t)  cnt                                                  term(\t or \n)
		return first_col_len + num_columns * (1 +       refresh::numeric_conversion_max_length<VALUE_T>()) + 1;
	}

public:
	MatrixOutputBuffer(TextFileWriter& writer, uint32_t first_col_len, size_t num_columns, size_t buff_size = 1ull << 23) :
		OutputBuffer(writer, get_max_record_len(first_col_len, num_columns), buff_size)
	{}

	void StoreKmer(const std::string& kmerSeq, const std::vector<VALUE_T>& values)
	{
		auto AsMatrixRow = [this](const std::string & kmerSeq, const std::vector<VALUE_T>&values, char* out) -> size_t
		{
			size_t res = store_kmer(kmerSeq, out, '\t');

			for (size_t i = 0; i < values.size() - 1; ++i)
				res += store_single_value(values[i], out, '\t');
			res += store_single_value(values.back(), out, '\n');

			return res;
		};

		shift_buffer_after_saving(AsMatrixRow(kmerSeq, values, get_buffer()));
	}

	void StoreKmer(const std::string& kmerSeq, const VALUE_T value)
	{
		auto AsMatrixRow_single_val = [this](const std::string & kmerSeq, const VALUE_T value, char* out) -> size_t
		{
			size_t res = store_kmer(kmerSeq, out, '\t');

			res += store_single_value(value, out, '\n');

			return res;
		};

		shift_buffer_after_saving(AsMatrixRow_single_val(kmerSeq, value, get_buffer()));
	}
};

class FastaOutputBuffer : public OutputBuffer
{
private:
	size_t get_max_record_len(uint32_t kmer_len) const
	{
		//     >\n k-mer                           \n
		return 2 + static_cast<size_t>(kmer_len) + 1;
	}

public:
	FastaOutputBuffer(TextFileWriter& writer, uint32_t kmer_len, size_t buff_size = 1ull << 23) :
		OutputBuffer(writer, get_max_record_len(kmer_len), buff_size)
	{}

	void StoreKmer(const std::string& kmerSeq)
	{
		auto AsFastaRecord = [this](const std::string & kmerSeq, char* out) -> size_t
		{
			out[0] = '>';
			out[1] = '\n';
			size_t res = 2;
			out += 2;

			res += store_kmer(kmerSeq, out, '\n');

			return res;
		};

		shift_buffer_after_saving(AsFastaRecord(kmerSeq, get_buffer()));
	}
};

inline void OutputBuffer::shift_buffer_after_saving(size_t record_len)
{
	if (out_buff_pos + max_line_len > buff_owner.size())
	{
		writer.Write(buff, out_buff_pos);
		out_buff_pos = 0;
	}

	out_buff_pos += record_len;
}

template<typename VALUE_T>
inline size_t OutputBuffer::store_single_value(const VALUE_T& val, char*& out, char term) const
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

inline size_t OutputBuffer::store_kmer(const std::string& kmerSeq, char*& out, char term) const
{
	std::memcpy(out, kmerSeq.data(), kmerSeq.length());
	out += kmerSeq.length();
	*out++ = term;

	return kmerSeq.length() + 1;
}
