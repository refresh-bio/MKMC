#include "FileGenerators.h"
#include <memory>



MatrixFileGenerator::MatrixFileGenerator(const Params& params, uint32_t binId, kmcdb::BinWriterSortedPlain<uint64_t>* bin) :
	bin(bin)
{
}



FASTAFileGenerator::FASTAFileGenerator(const Params& params, uint32_t binId, kmcdb::BinWriterSortedPlain<uint64_t>* /*bin*/) :
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
