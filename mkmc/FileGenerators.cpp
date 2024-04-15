#include "FileGenerators.h"
#include <memory>



MatrixFileGenerator::MatrixFileGenerator(const Params& params, uint32_t binId) :
	str_kmer_buff(std::make_unique<char[]>(params.mkmcParams.samples.size() * (params.mkmcParams.countSymbols + 1) + params.stage1Params.GetKmerLen() + 1)),
	file(params.mkmcParams.outputMatrixFiles[binId]),
	k(params.stage1Params.GetKmerLen())
{
	if (!file.is_open())
	{
		std::cerr << "Error: cannot create output file " << params.mkmcParams.outputMatrixFiles[binId] << "." << std::endl;
		exit(1);
	}

	file << "k-mer\t";
	for (const std::string& db : params.mkmcParams.samples)
	{
		file << db << '\t';
	}
	file << '\n';
}



FASTAFileGenerator::FASTAFileGenerator(const Params& params, uint32_t binId) :
	str_kmer_buff(std::make_unique<char[]>(params.stage1Params.GetKmerLen() + 1)),
	file(params.mkmcParams.outputFASTAFiles[binId]),
	k(params.stage1Params.GetKmerLen())
{
	if (!file.is_open())
	{
		std::cerr << "Error: cannot create output file " << params.mkmcParams.outputFASTAFiles[binId] << "." << std::endl;
		exit(1);
	}
}
