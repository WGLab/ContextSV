#ifndef SV_OBJECT_H
#define SV_OBJECT_H

#include <vector>
#include <memory>
#include <string>
#include <set>
#include <stdexcept>
#include <unordered_map>

#include "sv_types.h"

using namespace sv_types;

// Struct to represent a structural variant call
struct SVCall {
    uint32_t start;
    uint32_t end;
    // std::string sv_type = "NA";
    SVType sv_type = SVType::UNKNOWN;
    std::string alt_allele = ".";
    std::string data_type = "NA";
    std::string genotype = "./.";
    double hmm_likelihood = 0.0;
    int read_depth = 0;  // Breakpoint depth
    int support = 0;  // Number of supporting reads
    int cluster_size = 0;  // Number of SV calls in the cluster
    uint8_t qual = 0;  // Alignment quality score

    // Comparison operator for std::set
    bool operator<(const SVCall& other) const;

    // Constructor with parameters for all fields
    SVCall(uint32_t start, uint32_t end, SVType sv_type, const std::string& alt_allele, std::string data_type, std::string genotype, double hmm_likelihood, int read_depth, int support, int cluster_size, uint8_t qual) :
        start(start), end(end), sv_type(sv_type), alt_allele(alt_allele), data_type(data_type), genotype(genotype), hmm_likelihood(hmm_likelihood), read_depth(read_depth), support(support), cluster_size(cluster_size), qual(qual) {}
};

void addSVCall(std::vector<SVCall>& sv_calls, uint32_t start, uint32_t end, SVType sv_type, const std::string& alt_allele, std::string data_type, std::string genotype, double hmm_likelihood, int read_depth, uint8_t qual);

void mergeSVs(std::vector<SVCall>& sv_calls);

void filterSVsWithLowSupport(std::vector<SVCall> &sv_calls, int min_depth);

void filterSVsWithLowSupport(std::vector<SVCall> &sv_calls, int min_depth, const std::string& data_type);

uint32_t getSVCount(const std::vector<SVCall>& sv_calls);

void concatenateSVCalls(std::vector<SVCall>& sv_calls, const std::vector<SVCall>& sv_calls_update);

void mergeSVs(std::vector<SVCall> &sv_calls, double epsilon, int min_pts);

#endif // SV_OBJECT_H
