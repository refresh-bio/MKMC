#include "DumpWriter.h"



void DumpWriter::StoreHeader(const std::vector<std::string>& sample_names)
{
	fwrite("k-mer", 1, strlen("k-mer"), out);
	for (const auto& sample_name : sample_names)
	{
		fwrite("\t", 1, strlen("\t"), out);
		fwrite(sample_name.c_str(), 1, sample_name.length(), out);
	}
	fwrite("\n", 1, strlen("\n"), out);
}
