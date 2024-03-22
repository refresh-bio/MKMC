#include "FileGenerators.h"
#include <memory>



MatrixFileGenerator::MatrixFileGenerator(std::ostream& file, const std::vector<std::string>& samples, uint32_t countSymbols, uint32_t k) :
	str_kmer_buff(std::make_unique<char[]>(samples.size() * (countSymbols + 1) + k + 1)),
	file(file),
	k(k)
{
	file << "k-mer\t";
	for (const std::string& db : samples)
	{
		file << db << '\t';
	}
	file << '\n';
}



FASTAFileGenerator::FASTAFileGenerator(std::ostream& file, const std::vector<std::string>& samples, uint32_t countSymbols, uint32_t k) :
	str_kmer_buff(std::make_unique<char[]>(k + 1)),
	file(file),
	k(k)
{
}
