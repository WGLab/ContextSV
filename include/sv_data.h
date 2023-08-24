#ifndef SV_DATA_H
#define SV_DATA_H

#include "fasta_query.h"  // For querying the reference genome
#include "types.h"

#include <string>
#include <map>


class SVData {
    private:
        // SV candidate to read depth map
        SVDepthMap sv_calls;
        
        // Store a reference to the reference genome
        FASTAQuery *ref_genome;
        
    public:
        SVData(FASTAQuery& ref_genome);
        void addSVCall(std::string chr, int start, int end, int sv_type, std::string alt_allele);
        //void addSVCalls(SVData sv_calls);

        std::string getRefGenome();
        
        // Query the reference genome
        std::string getSequence(std::string chr, int pos_start, int pos_end);

        // Update the SV type for a given SV candidate
        void updateSVType(SVCandidate key, int sv_type);
        
        // Save SV calls to VCF
        void saveToVCF(FASTAQuery& ref_genome, std::string output_dir);

        // // Begin and end iterators for the SV candidate map
        // SVDepthMap::iterator begin() const { return this->sv_calls.begin(); }
        // SVDepthMap::iterator end() const { return this->sv_calls.end(); }

        // Begin and end iterators for the SV candidate map
        SVDepthMap::iterator begin() { return this->sv_calls.begin(); }
        SVDepthMap::iterator end() { return this->sv_calls.end(); }
};

#endif // SV_DATA_H
