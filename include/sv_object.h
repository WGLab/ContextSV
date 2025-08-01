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
    // SVDataType data_type = SVDataType::UNKNOWN;
    SVEvidenceFlags aln_type;
    Genotype genotype = Genotype::UNKNOWN;
    double hmm_likelihood = 0.0;
    int cn_state = 0;  // Copy number state
    int aln_offset = 0;  // Alignment offset (read vs. reference distance factor)
    int cluster_size = 0;  // Number of SV calls in the cluster

    bool operator<(const SVCall& other) const;

    SVCall() = default;

    SVCall(uint32_t start, uint32_t end, SVType sv_type, const std::string& alt_allele, SVEvidenceFlags aln_type, Genotype genotype, double hmm_likelihood, int cn_state, int aln_offset, int cluster_size) :
        start(start), end(end), sv_type(sv_type), alt_allele(alt_allele), aln_type(aln_type), genotype(genotype), hmm_likelihood(hmm_likelihood), cn_state(cn_state), aln_offset(aln_offset), cluster_size(cluster_size) {}
};

void addSVCall(std::vector<SVCall>& sv_calls, SVCall& sv_call);

// Merge SVs with identical start positions, and sum the cluster sizes
void mergeDuplicateSVs(std::vector<SVCall>& sv_calls);

uint32_t getSVCount(const std::vector<SVCall>& sv_calls);

// Merge SVs using DBSCAN clustering
void mergeSVs(std::vector<SVCall> &sv_calls, double epsilon, int min_pts, bool keep_noise, const std::string& json_filepath = "");

// Save clusters of SV calls to a JSON file
void saveClustersToJSON(const std::string& filename, const std::map<int, std::vector<SVCall>>& clusters);

#endif // SV_OBJECT_H
