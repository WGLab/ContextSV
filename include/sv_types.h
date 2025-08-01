#ifndef SV_TYPES_H
#define SV_TYPES_H

/// @cond
#include <string>
#include <map>
#include <set>
#include <unordered_map>
#include <tuple>
#include <bitset>
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
        LOH = 6  // Loss of heterozygosity
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
        {SVType::LOH, "LOH"}
    };

    // Mapping of SV types to symbols
    const std::unordered_map<SVType, std::string> SVTypeSymbol = {
        {SVType::UNKNOWN, "."},
        {SVType::DEL, "<DEL>"},
        {SVType::DUP, "<DUP>"},
        {SVType::INV, "<INV>"},
        {SVType::INS, "<INS>"},
        {SVType::BND, "<BND>"},
    };

    // Define constants for genotypes
    enum class Genotype {
        HOMOZYGOUS_REF = 0,
        HETEROZYGOUS = 1,
        HOMOZYGOUS_ALT = 2,
        UNKNOWN = 3
    };

    // Mapping of genotypes to strings
    const std::unordered_map<Genotype, std::string> GenotypeString = {
        {Genotype::HOMOZYGOUS_REF, "0/0"},
        {Genotype::HETEROZYGOUS, "0/1"},
        {Genotype::HOMOZYGOUS_ALT, "1/1"},
        {Genotype::UNKNOWN, "./."}
    };

    // Define constants for SV data types (evidence types)
    enum class SVDataType {
        CIGARINS = 0,
        CIGARDEL = 1,
        CIGARCLIP = 2,
        SPLIT = 3,
        SPLITDIST1 = 4,
        SPLITDIST2 = 5,
        SPLITINV = 6,
        SUPPINV = 7,
        HMM = 8,
        UNKNOWN = 9
    };

    using SVEvidenceFlags = std::bitset<10>;  // Bitset for SV data types

    // Mapping of SV data types to strings
    const std::unordered_map<SVDataType, std::string> SVDataTypeString = {
        {SVDataType::CIGARINS, "CIGARINS"},
        {SVDataType::CIGARDEL, "CIGARDEL"},
        {SVDataType::CIGARCLIP, "CIGARCLIP"},
        {SVDataType::SPLIT, "SPLIT"},
        {SVDataType::SPLITDIST1, "SPLITDIST1"},
        {SVDataType::SPLITDIST2, "SPLITDIST2"},
        {SVDataType::SPLITINV, "SPLITINV"},
        {SVDataType::SUPPINV, "SUPPINV"},
        {SVDataType::HMM, "HMM"},
        {SVDataType::UNKNOWN, "UNKNOWN"}
    };

    // Mapping of 6 copy number states to SV types
    const std::unordered_map<int, SVType> CNVTypeMap = {
        {0, SVType::UNKNOWN},
        {1, SVType::DEL},
        {2, SVType::DEL},
        {3, SVType::NEUTRAL},
        {4, SVType::LOH},
        {5, SVType::DUP},
        {6, SVType::DUP}
    };

    // Function to get the SV type string
    inline std::string getSVTypeString(SVType sv_type) {
        return SVTypeString.at(sv_type);
    }

    // Function to get the SV alignment type string from the bitset
    inline std::string getSVAlignmentTypeString(SVEvidenceFlags aln_type) {
        std::string result;
        for (size_t i = 0; i < SVDataTypeString.size(); ++i) {
            if (aln_type.test(i)) {
                result += SVDataTypeString.at(static_cast<SVDataType>(i)) + ",";
            }
        }
        if (!result.empty()) {
            result.pop_back();  // Remove the trailing comma
        }
        return result;
    }

    // Function to get the SV type from the CNV state
    inline SVType getSVTypeFromCNState(int cn_state) {
        return CNVTypeMap.at(cn_state);
    }

    // Function to get the genotype string
    inline std::string getGenotypeString(Genotype genotype) {
        return GenotypeString.at(genotype);
    }

    // Function to get the SV data type string
    // inline std::string getSVDataTypeString(SVDataType data_type) {
    //     return SVDataTypeString.at(data_type);
    // }

    // Function to get the SV type symbol
    inline std::string getSVTypeSymbol(SVType sv_type) {
        return SVTypeSymbol.at(sv_type);
    }

    // Function to check if an SV type is a valid update from copy number predictions
    inline bool isValidCopyNumberUpdate(SVType sv_type, SVType updated_sv_type) {
        if (updated_sv_type == SVType::UNKNOWN) {
            return false;
        } else if (sv_type == SVType::DEL && updated_sv_type != SVType::DEL) {
            return false;
        } else if (sv_type == SVType::INS && updated_sv_type != SVType::DUP) {
            return false;
        }
        return true;
    }
}

#endif // SV_TYPES_H
