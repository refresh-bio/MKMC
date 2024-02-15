#include "FileGenerators.h"
#include "../kmc/kmc_dump/nc_utils.h"
#include <memory>



MatrixFileGenerator::MatrixFileGenerator(std::ostream& file, const std::vector<std::string>& inputFiles, uint32_t countSymbols, uint32_t k) :
	str_kmer_buff(std::make_unique<char[]>(inputFiles.size() * (countSymbols + 1) + k + 1)),
	file(file),
	k(k)
{
	file << "k-mer\t";
	for (const std::string& db : inputFiles)
	{
		file << db << '\t';
	}
	file << '\n';
}

uint32_t MatrixFileGenerator::writeKmer(KMCFileWrapper::kmer_t& kmer, const std::vector<size_t>& kMersCounts)
{
	kmer.to_string(str_kmer_buff.get());
	uint32_t pos = k;
	for (size_t count : kMersCounts)
	{
		str_kmer_buff[pos++] = '\t';
		uint32_t shift = CNumericConversions::Int2PChar(count, reinterpret_cast<uchar*>(str_kmer_buff.get()) + pos);
		pos += shift;
	}

	str_kmer_buff[pos] = '\n';
	str_kmer_buff[pos + 1] = '\0';
	file << str_kmer_buff.get();

	return pos;
}

FASTAFileGenerator::FASTAFileGenerator(std::ostream& file, const std::vector<std::string>& inputFiles, uint32_t countSymbols, uint32_t k) :
	str_kmer_buff(std::make_unique<char[]>(k + 1)),
	file(file),
	k(k)
{
}

uint32_t FASTAFileGenerator::writeKmer(KMCFileWrapper::kmer_t& kmer, const std::vector<size_t>& kMersCounts)
{
	file << ">\n";
	kmer.to_string(str_kmer_buff.get());
	file << str_kmer_buff.get() << '\n';

	return k;
}
