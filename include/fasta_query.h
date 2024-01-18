// FASTAQuery: A class for querying a FASTA file.

#ifndef FASTA_QUERY_H
#define FASTA_QUERY_H

/// @cond
#include <string>
#include <map>
#include <unordered_map>
#include <vector>
/// @endcond

class FASTAQuery {
    private:
        std::string fasta_filepath;
        std::vector<std::string> chromosomes;
        std::unordered_map<std::string, std::string> chr_to_seq;

    public:
        int setFilepath(std::string fasta_filepath);
        std::string getFilepath();
        std::string query(std::string chr, int64_t pos_start, int64_t pos_end);

        // Get the chromosome contig lengths in VCF header format
        std::string getContigHeader();

        // Get the list of chromosomes, used for whole genome analysis
        std::vector<std::string> getChromosomes();

        // Check if a chromosome is in the FASTA file and return the chromosome name
        std::string hasChromosome(std::string chr);
};

#endif // FASTA_QUERY_H
