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
        SVData() {};

        int add(std::string chr, int64_t start, int64_t end, int sv_type, std::string alt_allele, std::string data_type, std::string genotype, double hmm_likelihood);

        void concatenate(const SVData& sv_data);

        // Update clipped base support for a given breakpoint location
        void updateClippedBaseSupport(std::string chr, int64_t pos);

        int getClippedBaseSupport(std::string chr, int64_t pos, int64_t end);
        
        void saveToVCF(FASTAQuery& ref_genome, std::string output_dir);

        std::map<SVCandidate, SVInfo>& getChromosomeSVs(std::string chr);

        std::set<std::string> getChromosomes();

        // Begin and end iterators for the SV candidate map
        SVDepthMap::iterator begin() { return this->sv_calls.begin(); }
        SVDepthMap::iterator end() { return this->sv_calls.end(); }

        // Get the total number of calls (For summary purposes)
        int totalCalls();
};

#endif // SV_DATA_H
