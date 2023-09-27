// FASTAQuery: A class for querying a FASTA file.

#ifndef FASTA_QUERY_H
#define FASTA_QUERY_H

/// @cond
#include <string>
#include <map>
/// @endcond

class FASTAQuery {
    private:
        std::string fasta_filepath;
        std::map<std::string, std::string> chr_to_seq;

    public:
        int setFilepath(std::string fasta_filepath);
        std::string getFilepath();
        std::string query(std::string chr, int pos_start, int pos_end);

        // Get the chromosome contig lengths in VCF header format
        std::string getContigHeader();
};

#endif // FASTA_QUERY_H
