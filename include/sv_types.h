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
    static const int UNKNOWN = -1;

    // Define SVTypeString for SV types
    static const std::string SVTypeString[] = {"DEL", "DUP", "INV", "INS", "BND", "UNKNOWN"};

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

    // SV (start, end, alt_allele)
    using SVCandidate = std::tuple<int64_t, int64_t, std::string>;
    
    // Chromosome to SV candidate to read depth map
    using SVDepthMap = std::unordered_map<std::string, std::map<SVCandidate, SVInfo>>;

    // Define a map for storing copy number calls by SV candidate
    using SVCopyNumberMap = std::map<SVCandidate, std::tuple<int, std::string, std::string>>;
}

#endif // SV_TYPES_H
