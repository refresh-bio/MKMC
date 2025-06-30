#pragma once
#include <string>

//#include "refresh/archive/lib/archive_output.h" - has to be included by kmcdb
//#include "refresh/archive/lib/archive_input.h" - has to be included by kmcdb
#include "kmcdb/kmcdb.h"

class MatrixStatsWriter
{
	refresh::archive_output archive;
	bool is_closed = false;
	public:
	MatrixStatsWriter(const std::string& path)
	{
		if (!archive.open_file_buffered(path, false))
			throw std::runtime_error("Cannot open output archive: " + path);
	}
	//each add must be of unique name
	void Add(const std::string& name, const std::vector<uint8_t>& data)
	{
		assert(!is_closed);
		if (archive.get_stream_id(name) != -1)
			throw std::runtime_error("Each stream name in MatrixStatsWriter::Add must be unique, but this is not: " + name);
		auto id = archive.register_stream(name);
		if (id == -1)
			throw std::runtime_error("Cannot register stream: " + name);

		if (!archive.add_part(id, data))
			throw std::runtime_error("Cannot store data in stream: " + name);
	}

	void Close()
	{
		assert(!is_closed);
		if(!archive.close())
			throw std::runtime_error("Cannot close archive in store data in MatrixStatsWriter");
		is_closed = true;
	}
	//Scott Meyers, Effective C++, Item 11: Prevent exceptions from leaving destructors.
	~MatrixStatsWriter() noexcept
	{
		if (is_closed)
			return;
		try
		{
			Close();
		}
		catch (...) {}
	}
};

class MatrixStatsReader
{
	refresh::archive_input archive;
	bool is_closed = false;
public:
	MatrixStatsReader(const std::string& path)
	{
		if (!archive.open_file_buffered(path, false))
			throw std::runtime_error("Cannot open output archive: " + path);
	}
	bool HasStream(const std::string& name)
	{
		assert(!is_closed);
		return archive.get_stream_id(name) != -1;
	}
	bool Get(const std::string& name, std::vector<uint8_t>& data)
	{
		assert(!is_closed);
		data.clear();
		auto id = archive.get_stream_id(name);
		if (id == -1)
			return false;
		uint64_t meta;
		if (!archive.get_part(id, data, meta))
			throw std::runtime_error("Cannot read data from stream: " + name);

		return true;
	}
	void Close()
	{
		assert(!is_closed);
		if (!archive.close())
			throw std::runtime_error("Cannot close archive in store data in MatrixStatsWriter");
		is_closed = true;
	}
	//Scott Meyers, Effective C++, Item 11: Prevent exceptions from leaving destructors.
	~MatrixStatsReader() noexcept
	{
		if (is_closed)
			return;
		try
		{
			Close();
		}
		catch (...) {}
	}
};