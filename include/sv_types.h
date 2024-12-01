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
}

#endif // SV_TYPES_H
