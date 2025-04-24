#include "TextFileWritingUtilities.h"



void TextFileWriter::StoreHeader(const std::vector<std::string>& data_columns_headers, const std::string& first_column_header)
{
	fwrite(first_column_header.c_str(), 1, first_column_header.length(), out);
	for (const auto& column_name : data_columns_headers)
	{
		fwrite("\t", 1, strlen("\t"), out);
		fwrite(column_name.c_str(), 1, column_name.length(), out);
	}
	fwrite("\n", 1, strlen("\n"), out);
}

OutputBuffer::~OutputBuffer()
{
	if (out_buff_pos)
	{
		writer.Write(buff, out_buff_pos);
		out_buff_pos = 0;
	}
}
