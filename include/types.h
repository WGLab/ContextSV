#ifndef TYPES_H
#define TYPES_H

#include <string>
#include <vector>
#include <map>

// CNV candidate location map
// (chr, snp_pos) : cnv_type
using SNPLocation = std::pair<std::string, int>;
using SNPToCNVMap = std::map<SNPLocation, int>;

// SV candidate location map
// (chr, start, end) : (sv_type, read_depth, ref_allele, alt_allele)
using SVCandidate = std::tuple<std::string, int, int>;
using SVInfo = std::tuple<int, int, std::string, std::string>;
using SVInfoMap = std::map<SVCandidate, SVInfo>;

#endif // TYPES_H
