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
    enum class SVType {
        UNKNOWN = -1,
        DEL = 0,
        DUP = 1,
        INV = 2,
        INS = 3,
        BND = 4,
        NEUTRAL = 5,  // Neutral copy number with unknown type
        INV_DUP = 6,  // Inversion duplication
        COMPLEX = 7  // Complex SV
    };

    // Mapping of SV types to strings
    const std::unordered_map<SVType, std::string> SVTypeString = {
        {SVType::UNKNOWN, "UNKNOWN"},
        {SVType::DEL, "DEL"},
        {SVType::DUP, "DUP"},
        {SVType::INV, "INV"},
        {SVType::INS, "INS"},
        {SVType::BND, "BND"},
        {SVType::NEUTRAL, "NEUTRAL"},
        {SVType::INV_DUP, "INVDUP"},
        {SVType::COMPLEX, "COMPLEX"}
    };

    // Mapping of 6 copy number states to SV types
    const std::unordered_map<int, SVType> CNVTypeMap = {
        {0, SVType::UNKNOWN},
        {1, SVType::DEL},
        {2, SVType::DEL},
        {3, SVType::NEUTRAL},
        {4, SVType::NEUTRAL},
        {5, SVType::DUP},
        {6, SVType::DUP}
    };

    // Function to get the SV type string
    inline std::string getSVTypeString(SVType sv_type) {
        return SVTypeString.at(sv_type);
    }

    // Function to get the SV type from the CNV state
    inline SVType getSVTypeFromCNState(int cn_state) {
        return CNVTypeMap.at(cn_state);
    }

    // Create a struct for storing SV information
    struct SVInfo {
        SVType sv_type;
        int read_support;  // Number of reads supporting the SV breakpoints
        int read_depth;  // Read depth at the SV start position
        std::set<std::string> data_type;  // Alignment type used to call the SV
        int sv_length;
        std::string genotype = "./.";  // Default genotype (no call)
        double hmm_likelihood = 0.0;  // HMM likelihood score for the state sequence

        SVInfo() = default;
        SVInfo(SVType sv_type, int read_support, int read_depth, std::string data_type, int sv_length, std::string genotype, double hmm_likelihood) :
            sv_type(sv_type), read_support(read_support), read_depth(read_depth), data_type({data_type}), sv_length(sv_length), genotype(genotype), hmm_likelihood(hmm_likelihood) {}
    };

    // Type definition for SV-related structures
    using SVCandidate = std::tuple<int32_t, int32_t, std::string>;  // SV (start, end, alt_allele)
    using SVDepthMap = std::unordered_map<std::string, std::map<SVCandidate, SVInfo>>;  // Chromosome -> SV candidate -> SV info
}

#endif // SV_TYPES_H
