// This file contains all the type definitions used in the project

#ifndef TYPES_H
#define TYPES_H

/// @cond
#include <string>
#include <vector>
#include <map>
#include <set>
/// @endcond

// SNP data is a struct containing vectors used in predicting copy number states
struct SNPData {
    std::vector<int64_t> locations;
    std::vector<double> pfbs;
    std::vector<double> bafs;
    std::vector<double> log2_ratios;
    std::vector<int> state_sequence;

    SNPData():
        locations({}),\
        pfbs({}), \
        bafs({}), \
        log2_ratios({}), \
        state_sequence({}) {}
};

// CNV candidate location map
// (chr, snp_pos) : cnv_type
using SNPLocation = std::pair<std::string, int64_t>;
using SNPToCNVMap = std::map<SNPLocation,  int>;

// SV candidate read depth map. An SV is defined by its location, type, and
// alternate allele.
// (chr, start, end, sv_type, alt_allele) : (ref_allele, num_reads)
using SVCandidate = std::tuple<std::string, int64_t, int64_t, std::string>;  // chr, start, end, alt_allele

// Create a struct for storing SV information
struct SVInfo {
    int sv_type;
    int read_depth;
    std::set<std::string> data_type;
    int sv_length;

    SVInfo() :
        sv_type(-1), read_depth(0), data_type({}), sv_length(0) {}
        
    SVInfo(int sv_type, int read_depth, std::string data_type, int sv_length) :
        sv_type(sv_type), read_depth(read_depth), data_type({data_type}), sv_length(sv_length) {}
};

using SVDepthMap = std::map<SVCandidate, SVInfo>;  // Map for getting type and read depth of SV candidates

// Alignment data used for split-read SV calling (chr, start, end, sequence)
using AlignmentData   = std::tuple<std::string, int64_t, int64_t, std::string>;
using AlignmentVector = std::vector<AlignmentData>;

// Query map (query name, alignment vector)
using QueryMap = std::map<std::string, AlignmentVector>;

#endif // TYPES_H
