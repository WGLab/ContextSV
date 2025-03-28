#ifndef SV_OBJECT_H
#define SV_OBJECT_H

#include <vector>
#include <memory>
#include <string>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <string_view>

#include "sv_types.h"

using namespace sv_types;

struct SVCall {
    uint32_t start = 0;
    uint32_t end = 0;
    SVType sv_type = SVType::UNKNOWN;
    std::string alt_allele = ".";
    SVDataType data_type = SVDataType::UNKNOWN;
    Genotype genotype = Genotype::UNKNOWN;
    double hmm_likelihood = 0.0;
    int read_depth = 0;  // Breakpoint depth
    double mismatch_rate = 0.0;  // Highest mismatch rate in reads used for the SV call
    int cluster_size = 0;  // Number of SV calls in the cluster

    bool operator<(const SVCall& other) const;

    SVCall() = default;

    SVCall(uint32_t start, uint32_t end, SVType sv_type, std::string alt_allele, SVDataType data_type, Genotype genotype, double hmm_likelihood, int read_depth, double mismatch_rate, int cluster_size) :
        start(start), end(end), sv_type(sv_type), alt_allele(alt_allele), data_type(data_type), genotype(genotype), hmm_likelihood(hmm_likelihood), read_depth(read_depth), mismatch_rate(mismatch_rate), cluster_size(cluster_size) {}
};

void addSVCall(std::vector<SVCall>& sv_calls, SVCall& sv_call);

// Merge SVs with identical start positions, and sum the cluster sizes
void mergeDuplicateSVs(std::vector<SVCall>& sv_calls);

uint32_t getSVCount(const std::vector<SVCall>& sv_calls);

void concatenateSVCalls(std::vector<SVCall>& sv_calls, const std::vector<SVCall>& sv_calls_update);

// Merge SVs using DBSCAN clustering
void mergeSVs(std::vector<SVCall> &sv_calls, double epsilon, int min_pts, bool keep_noise);

#endif // SV_OBJECT_H
