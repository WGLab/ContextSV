#ifndef SV_DATA_H
#define SV_DATA_H


/// @cond
#include <string>
#include <map>
#include <set>
#include <mutex>

#include "sv_types.h"
#include "fasta_query.h"
/// @endcond

// Include the SV types namespace
using namespace sv_types;

// SV data class
class SVData {
    private:
        SVDepthMap sv_calls;

        // Map of clipped base support by position (chr, pos) : depth
        std::map<std::pair<std::string, int64_t>, int> clipped_base_support;
        
    public:
        SVData() {};

        int add(std::string chr, int32_t start, int32_t end, SVType sv_type, std::string alt_allele, std::string data_type, std::string genotype, double hmm_likelihood);

        void concatenate(const SVData& sv_data);

        // Update clipped base support for a given breakpoint location
        void updateClippedBaseSupport(std::string chr, int64_t pos);

        int getClippedBaseSupport(std::string chr, int64_t pos, int64_t end);
        
        void saveToVCF(ReferenceGenome& ref_genome, std::string output_dir);

        std::map<SVCandidate, SVInfo>& getChromosomeSVs(std::string chr);

        std::set<std::string> getChromosomes();

        // Begin and end iterators for the SV candidate map
        SVDepthMap::iterator begin() { return this->sv_calls.begin(); }
        SVDepthMap::iterator end() { return this->sv_calls.end(); }

        // Get the total number of calls (For summary purposes)
        int totalCalls();
};

#endif // SV_DATA_H
