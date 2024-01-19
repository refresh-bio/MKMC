#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cinttypes>
#include <limits>
#include <memory>
#include <sstream>
#include "../kmc/kmc_api/kmer_api.h"
#include "../kmc/kmc_api/kmc_file.h"
#include "Dump.h"
#include "Filter.h"



bool Dump::allKAreSame(const std::vector<KMCFileWrapper>& samples)
{
    if (samples.empty())
        return true;
    uint32_t k = samples.front().GetK();
    for (const auto& sample : samples)
        if (k != sample.GetK())
            return false;
    return true;
}

void Dump::openDatabases(std::vector<KMCFileWrapper>& samples)
{
    
    for (const std::string& fileName : params.mkmcParams.toolsOutputFiles)
    {
        samples.emplace_back(fileName);
    }

    if (!allKAreSame(samples))
    {
        std::cerr << "Error: each database should have the same k\n";
        exit(1);
    }
}

void Dump::dumpToStd()
{
    std::vector<KMCFileWrapper> samples;
    openDatabases(samples);

    std::string str_kmer;
    uint64_t dumped_symbols = 0;
    std::vector<size_t> kMersCounts(samples.size());
    Filter filter(params);
    while (true)
    {
        std::size_t min_id = std::numeric_limits<size_t>::max();
        for (std::size_t i = 0; i < samples.size(); ++i)
        {
            if (!samples[i].Finished())
            {
                if (min_id == std::numeric_limits<size_t>::max() || samples[i].First() < samples[min_id].First())
                    min_id = i;
            }
        }
        if (min_id == std::numeric_limits<size_t>::max()) // no more k-mers
            break;

        auto min_kmer = samples[min_id].First();
        min_kmer.to_string(str_kmer);

        for (size_t itSample = 0; itSample < samples.size(); ++itSample)
        {
            KMCFileWrapper& sample = samples[itSample];
            size_t count;
            if (sample.Finished() || !(sample.First() == min_kmer))
                count = 0;
            else
            {
                count = sample.FirstCount();
                sample.Next();
            }

            kMersCounts[itSample] = count;
        }

        if (filter.keepKMer(kMersCounts))
        {
            std::cout << str_kmer;
            dumped_symbols += str_kmer.length();

            std::ostringstream sstream;
            for (size_t count : kMersCounts)
            {
                sstream << "\t" << count;
            }

            std::cout << sstream.str();
            std::cout << "\n";
            dumped_symbols += sstream.str().length() + 1;

            if (dumped_symbols >= params.mkmcParams.dumpStepSize)
            {
                std::system("pause");
                dumped_symbols = 0;
            }
        }
    }
}

void Dump::dumpToFile()
{
    std::ofstream output_file(params.mkmcParams.outputFile);
    if (!output_file.is_open())
    {
        std::cerr << "Error: cannot create output file " << params.mkmcParams.outputFile << "\n";
        exit(1);
    }

    std::vector<KMCFileWrapper> samples;
    openDatabases(samples);

    std::string str_kmer;
    std::vector<size_t> kMersCounts(samples.size());
    Filter filter(params);
    while (true)
    {
        std::size_t min_id = std::numeric_limits<size_t>::max();
        for (std::size_t i = 0; i < samples.size(); ++i)
        {
            if (!samples[i].Finished())
            {
                if (min_id == std::numeric_limits<size_t>::max() || samples[i].First() < samples[min_id].First())
                    min_id = i;
            }
        }
        if (min_id == std::numeric_limits<size_t>::max()) // no more k-mers
            break;

        auto min_kmer = samples[min_id].First();
        min_kmer.to_string(str_kmer);

        for (size_t itSample = 0; itSample < samples.size(); ++itSample)
        {
            KMCFileWrapper& sample = samples[itSample];
            size_t count;
            if (sample.Finished() || !(sample.First() == min_kmer))
                count = 0;
            else
            {
                count = sample.FirstCount();
                sample.Next();
            }

            kMersCounts[itSample] = count;
        }

        if (filter.keepKMer(kMersCounts))
        {
            output_file << str_kmer;

            std::ostringstream sstream;
            for (size_t count : kMersCounts)
            {
                output_file << "\t" << count;
            }

            output_file << "\n";
        }
    }
}

KMCFileWrapper::KMCFileWrapper(const std::string& path)
{
    kmc_file = std::make_unique<CKMCFile>();
    if (!kmc_file->OpenForListing(path))
    {
        std::cerr << "Error: cannot open kmc database " << path << "\n";
        exit(1);
    }
    if (kmc_file->IsKMC2() == true)
    {
        std::cerr << "Error: kmc database not sorted: " << path << "\n";
        std::cerr << "Hint: KMC database may be sorted with kmc_tools transform <input> sort <output>\n";
        exit(1);
    }
    CKMCFileInfo kmc_file_info;
    kmc_file->Info(kmc_file_info);
    k = kmc_file_info.kmer_length;
    tot_kmers = kmc_file_info.total_kmers;

    cur = kmer_t(k);

    if (!Finished())
        Next();
}

inline void KMCFileWrapper::Next()
{
    if (!kmc_file->ReadNextKmer(cur, cur_count))
    {
        // for the last one k-mer (cur_kmer_no == tot_kmers), it is legal to call Next, but ReadNextKmer will fail
        if (cur_kmer_no != tot_kmers)
        {
            std::cerr << "Error: critical, this should not happen, details: " << __FILE__ << "(" << __LINE__ << ")\n";
            exit(1);
        }
    }
    ++cur_kmer_no;
}
