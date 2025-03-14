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

struct SVCall {
    uint32_t start;
    uint32_t end;
    SVType sv_type = SVType::UNKNOWN;
    std::string alt_allele = ".";
    std::string data_type = "NA";
    std::string genotype = "./.";
    double hmm_likelihood = 0.0;
    int read_depth = 0;  // Breakpoint depth
    int support = 0;  // Number of supporting reads
    int cluster_size = 0;  // Number of SV calls in the cluster

    bool operator<(const SVCall& other) const;

    SVCall(uint32_t start, uint32_t end, SVType sv_type, const std::string& alt_allele, std::string data_type, std::string genotype, double hmm_likelihood, int read_depth, int support, int cluster_size) :
        start(start), end(end), sv_type(sv_type), alt_allele(alt_allele), data_type(data_type), genotype(genotype), hmm_likelihood(hmm_likelihood), read_depth(read_depth), support(support), cluster_size(cluster_size) {}
};

void addSVCall(std::vector<SVCall>& sv_calls, SVCall& sv_call);

// Merge SVs with identical start positions, and sum the cluster sizes
void mergeDuplicateSVs(std::vector<SVCall>& sv_calls);

uint32_t getSVCount(const std::vector<SVCall>& sv_calls);

void concatenateSVCalls(std::vector<SVCall>& sv_calls, const std::vector<SVCall>& sv_calls_update);

// Merge SVs using DBSCAN clustering
void mergeSVs(std::vector<SVCall> &sv_calls, double epsilon, int min_pts);

#endif // SV_OBJECT_H
