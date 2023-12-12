#include <iostream>
#include <fstream>
#include <string>
#include  "kmc_core/kmc_runner.h"
#include  "kmc_api/kmc_file.h"
#include  "kmc_api/kmer_api.h"

void generate_test_file(const std::string& path)
{
	std::string seq = "AGCTACTACTGACTGACTTACTATGCTGATCGTACACACATGAC";
	std::ofstream out(path, std::ios::binary);
	if (!out)
	{
		std::cerr << "Error: cannot open file" << path << "\n";
		exit(1);
	}
	out.write(">\n", 2);
	out.write(seq.c_str(), seq.length());
	out.write("\n", 1);

}

int main(int argc, char**argv)
{
    using namespace std::string_literals;

    auto path = "test.fa"s;
	generate_test_file(path);
	int k = 20;
    auto db_path = std::to_string(k) + "-mers";
    try
    {
        KMC::Runner runner;

        KMC::Stage1Params stage1Params;
        stage1Params
            .SetKmerLen(k)
			.SetInputFileType(KMC::InputFileType::FASTA)
            .SetInputFiles({ path });

        auto stage1Result = runner.RunStage1(stage1Params);

        KMC::Stage2Params stage2Params;

        stage2Params
			.SetCutoffMin(1)
            .SetOutputFileName(db_path);

        auto stage2Result = runner.RunStage2(stage2Params);

        //print some stats
        std::cout << "total k-mers: " << stage2Result.nTotalKmers << "\n";
        std::cout << "total unique k-mers: " << stage2Result.nUniqueKmers << "\n";
    }
    catch (const std::exception& e)
    {
        std::cerr << e.what() << '\n';
    }

    CKMCFile kmc_file;
    if(!kmc_file.OpenForListing(db_path))
    {
        std::cerr << "Error: cannot open kmc database: " << db_path << "\n";
        return  EXIT_FAILURE;
    }
 
    CKmerAPI kmer(k);
    uint64 count;
    while(kmc_file.ReadNextKmer(kmer, count))
        std::cerr << kmer.to_string() << "\t" << count << "\n";
}