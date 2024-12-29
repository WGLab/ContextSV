// ReferenceGenome: A class for querying a reference genome FASTA file.

#ifndef FASTA_QUERY_H
#define FASTA_QUERY_H

/// @cond
#include <string>
#include <map>
#include <unordered_map>
#include <vector>
#include <mutex>
/// @endcond

class ReferenceGenome {
    private:
        std::string fasta_filepath;
        std::vector<std::string> chromosomes;
        std::unordered_map<std::string, std::string> chr_to_seq;
        //mutable std::mutex mtx;
        std::mutex& shared_mutex;

    public:
	    ReferenceGenome(std::mutex& mtx) : shared_mutex(mtx) {}
    
        int setFilepath(std::string fasta_filepath);
        std::string getFilepath() const;
        std::string query(const std::string& chr, uint32_t pos_start, uint32_t pos_end) const;

        // Get the chromosome contig lengths in VCF header format
        std::string getContigHeader() const;

        // Get the list of chromosomes, used for whole genome analysis
        std::vector<std::string> getChromosomes() const;

        // Get the length of a chromosome
        uint32_t getChromosomeLength(std::string chr) const;
};

#endif // FASTA_QUERY_H
