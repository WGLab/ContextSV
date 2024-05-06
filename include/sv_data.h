#ifndef SV_DATA_H
#define SV_DATA_H

#include "fasta_query.h"  // For querying the reference genome

/// @cond
#include <string>
#include <map>
#include <set>
#include <mutex>

#include "sv_types.h"
/// @endcond

// Include the SV types namespace
using namespace sv_types;



// SV data class
class SVData {
    private:
        SVDepthMap sv_calls;

        // Map of clipped base support by position (chr, pos) : depth
        std::map<std::pair<std::string, int64_t>, int> clipped_base_support;

        // SV type to string map for VCF output
        std::map<int, std::string> sv_type_map = {
            {0, "DEL"},
            {1, "DUP"},
            {2, "INV"},
            {3, "INS"},
            {4, "BND"},
            {5, "DUP"}
        };
        
    public:
        // Constructor
        SVData() {};

        // Add a new SV candidate to the map
        int add(std::string chr, int64_t start, int64_t end, int sv_type, std::string alt_allele, std::string data_type);

        // Concatenates two SVData objects
        void concatenate(const SVData& sv_data);

        // Update the SV type for a given SV candidate
        void updateSVType(std::string chr, SVCandidate key, int sv_type, std::string data_type);

        // Update the SV genotype for a given SV candidate
        void updateGenotype(std::string chr, SVCandidate key, std::string genotype);

        // Update clipped base support for a given breakpoint location
        void updateClippedBaseSupport(std::string chr, int64_t pos);

        // Get the SV clipped base support
        int getClippedBaseSupport(std::string chr, int64_t pos, int64_t end);
        
        // Save SV calls to VCF
        void saveToVCF(FASTAQuery& ref_genome, std::string output_dir);

        // Get the chromosome SVs
        std::map<SVCandidate, SVInfo>& getChromosomeSVs(std::string chr);

        // Get the chromosomes
        std::set<std::string> getChromosomes();

        // Begin and end iterators for the SV candidate map
        SVDepthMap::iterator begin() { return this->sv_calls.begin(); }
        SVDepthMap::iterator end() { return this->sv_calls.end(); }

        // Get the total number of calls (For summary purposes)
        int totalCalls();

        // Get the total number of deletions (For testing purposes)
        int totalDeletions();

        // Get the total number of deletions for a chromosome (For testing purposes)
        int totalDeletions(std::string chr);
};

#endif // SV_DATA_H
