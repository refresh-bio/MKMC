#include "DumpWriter.h"



void DumpWriter::StoreHeader(const std::vector<std::string>& sample_names, const std::string& first_item)
{
	fwrite(first_item.c_str(), 1, first_item.length(), out);
	for (const auto& sample_name : sample_names)
	{
		fwrite("\t", 1, strlen("\t"), out);
		fwrite(sample_name.c_str(), 1, sample_name.length(), out);
	}
	fwrite("\n", 1, strlen("\n"), out);
}
