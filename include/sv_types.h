#ifndef SV_TYPES_H
#define SV_TYPES_H

/// @cond
#include <string>
#include <map>
#include <set>
#include <unordered_map>
#include <tuple>
/// @endcond

namespace sv_types {
    // Define constants for SV types
    static const int DEL = 0;
    static const int DUP = 1;
    static const int INV = 2;
    static const int INS = 3;
    static const int BND = 4;
    static const int TANDUP = 5;  // Tandem duplication
    static const int UNKNOWN = -1;

    // Define SVTypeString for SV types
    static const std::string SVTypeString[] = {"DEL", "DUP", "INV", "INS", "BND", "DUP"};

    // Create a struct for storing SV information
    struct SVInfo {
        int sv_type;
        int read_support;  // Number of reads supporting the SV breakpoints
        int read_depth;  // Read depth at the SV start position
        std::set<std::string> data_type;  // Alignment type used to call the SV
        int sv_length;
        std::string genotype = "./.";  // Default genotype (no call)
        double hmm_likelihood = 0.0;  // HMM likelihood score for the state sequence
        std::string read_id = "";  // Read ID supporting the SV for debugging (just keep the first one found)

        SVInfo() :
            sv_type(-1), read_support(0), read_depth(0), data_type({}), sv_length(0), genotype("./."), hmm_likelihood(0.0), read_id("") {}
            
        SVInfo(int sv_type, int read_support, int read_depth, std::string data_type, int sv_length, std::string genotype, double hmm_likelihood, std::string read_id) :
            sv_type(sv_type), read_support(read_support), read_depth(read_depth), data_type({data_type}), sv_length(sv_length), genotype(genotype), hmm_likelihood(hmm_likelihood), read_id(read_id) {}
    };

    // SV (start, end, alt_allele)
    using SVCandidate = std::tuple<int64_t, int64_t, std::string>;
    
    // Chromosome to SV candidate to read depth map
    using SVDepthMap = std::unordered_map<std::string, std::map<SVCandidate, SVInfo>>;

    // Define a map for storing copy number calls by SV candidate
    using SVCopyNumberMap = std::map<SVCandidate, std::tuple<int, std::string, std::string>>;

    // Create a type for storing SV update information from copy number caller
    // (SVCandidate, SV type, genotype, data type)
    using SVUpdate = std::tuple<SVCandidate, int, std::string, std::string>;
}

#endif // SV_TYPES_H
