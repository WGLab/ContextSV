// This file contains all the type definitions used in the project

#ifndef TYPES_H
#define TYPES_H

#include <string>
#include <vector>
#include <map>

// CNV candidate location map
// (chr, snp_pos) : cnv_type
using SNPLocation = std::pair<std::string, int>;
using SNPToCNVMap = std::map<SNPLocation, int>;

// SV candidate read depth map. An SV is defined by its location, type, and
// alternate allele.
// (chr, start, end, sv_type, alt_allele) : (ref_allele, num_reads)
using SVCandidate = std::tuple<std::string, int, int, std::string>;  // chr, start, end, alt_allele
using SVInfo = std::pair<int, int>;  // SV type, read depth
using SVDepthMap = std::map<SVCandidate, SVInfo>;  // Map for getting type and read depth of SV candidates

// SV calling:
// Alignment location (chr, start, end, depth)
using AlignmentData = std::tuple<std::string, int, int, int>;
using AlignmentVector = std::vector<AlignmentData>;

// Query map (query name, alignment vector)
using QueryMap = std::map<std::string, AlignmentVector>;

#endif // TYPES_H
