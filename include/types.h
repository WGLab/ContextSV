// This file contains all the type definitions used in the project

#ifndef TYPES_H
#define TYPES_H

/// @cond
#include <string>
#include <vector>
#include <map>
#include <set>
/// @endcond

// CNV candidate location map
// (chr, snp_pos) : cnv_type
using SNPLocation = std::pair<std::string, int>;
using SNPToCNVMap = std::map<SNPLocation, int>;

// SV candidate read depth map. An SV is defined by its location, type, and
// alternate allele.
// (chr, start, end, sv_type, alt_allele) : (ref_allele, num_reads)
using SVCandidate = std::tuple<std::string, int, int, std::string>;  // chr, start, end, alt_allele
// using SVInfo = std::tuple<int, int, std::string>;  // SV type, read depth, alignment type (CIGAR or split read)

// // Create a class for storing SV information
// class SVInfo {
//     public:
//         int sv_type;
//         int read_depth;
//         std::set<std::string> data_type;

//         SVInfo(int sv_type, int read_depth, std::string data_type) : sv_type(sv_type), read_depth(read_depth), data_type({data_type}) {}

//         SVInfo() : sv_type(-1), read_depth(0) {}
// };
// Create a struct for storing SV information
struct SVInfo {
    int sv_type;
    int read_depth;
    std::set<std::string> data_type;

    SVInfo() : sv_type(-1), read_depth(0) {}
    SVInfo(int sv_type, int read_depth, std::string data_type) : sv_type(sv_type), read_depth(read_depth), data_type({data_type}) {}
};

using SVDepthMap = std::map<SVCandidate, SVInfo>;  // Map for getting type and read depth of SV candidates

// SV calling:
// Alignment location (chr, start, end, depth)
using AlignmentData = std::tuple<std::string, int, int, int>;
using AlignmentVector = std::vector<AlignmentData>;

// Query map (query name, alignment vector)
using QueryMap = std::map<std::string, AlignmentVector>;

#endif // TYPES_H
