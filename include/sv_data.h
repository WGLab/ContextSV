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

// Create a struct for storing SV information
struct SVInfo {
    int sv_type;
    int read_depth;
    std::set<std::string> data_type;  // Alignment type used to call the SV
    int sv_length;
    std::string genotype = "./.";  // Default genotype (no call)

    SVInfo() :
        sv_type(-1), read_depth(0), data_type({}), sv_length(0), genotype("./.") {}
        
    SVInfo(int sv_type, int read_depth, std::string data_type, int sv_length, std::string genotype) :
        sv_type(sv_type), read_depth(read_depth), data_type({data_type}), sv_length(sv_length), genotype(genotype) {}
};

using SVCandidate = std::tuple<std::string, int64_t, int64_t, std::string>;  // SV (chr, start, end, alt_allele)
using SVDepthMap = std::map<SVCandidate, SVInfo>;  // SV candidate to read depth map

// SV data class
class SVData {
    private:
        SVDepthMap sv_calls;
        // mutable std::mutex sv_calls_mtx;  // Mutex for locking the SV candidate map

        // Map of clipped base support by position (chr, pos) : depth
        std::map<std::pair<std::string, int64_t>, int> clipped_base_support;

        // SV type to string map (DEL, INS, INV, DUP, BND)
        // DUPs [1] are INS with INFO/REPTYPE=DUP
        std::map<int, std::string> sv_type_map = {
            {0, "DEL"},
            {1, "DUP"},
            {2, "INV"},
            {3, "INS"},
            {4, "BND"}
        };
        
    public:
        // // Define constants for SV types
        // static const int DEL = 0;
        // static const int DUP = 1;
        // static const int INV = 2;
        // static const int INS = 3;
        // static const int BND = 4;
        // static const int UNKNOWN = -1;

        // Constructor
        SVData() {};

        // Add a new SV candidate to the map
        void add(std::string chr, int64_t start, int64_t end, int sv_type, std::string alt_allele, std::string data_type);

        // Concatenates two SVData objects
        void concatenate(const SVData& sv_data);

        // Update the SV type for a given SV candidate
        void updateSVType(SVCandidate key, int sv_type, std::string data_type);

        // Update the SV genotype for a given SV candidate
        void updateGenotype(SVCandidate key, std::string genotype);

        // Update clipped base support for a given breakpoint location
        void updateClippedBaseSupport(std::string chr, int64_t pos);

        // Get the SV clipped base support
        int getClippedBaseSupport(std::string chr, int64_t pos, int64_t end);
        
        // Save SV calls to VCF
        void saveToVCF(FASTAQuery& ref_genome, std::string output_dir);

        // Begin and end iterators for the SV candidate map
        SVDepthMap::iterator begin() { return this->sv_calls.begin(); }
        SVDepthMap::iterator end() { return this->sv_calls.end(); }

        // Define the size of the SV candidate map
        int size() { return this->sv_calls.size(); }
};

#endif // SV_DATA_H
