#ifndef SV_OBJECT_H
#define SV_OBJECT_H

#include <vector>
#include <memory>
#include <string>
#include <set>
#include <stdexcept>
#include <unordered_map>

// Struct to represent a structural variant call
struct SVCall {
    uint32_t start;
    uint32_t end;
    std::string sv_type = "NA";
    std::string alt_allele = ".";
    std::string data_type = "NA";
    std::string genotype = "./.";
    double hmm_likelihood = 0.0;
    int support = 0;  // Exact breakpoint support
    int total_support = 0;  // Support at either breakpoint

    // Comparison operator for std::set
    bool operator<(const SVCall& other) const;

    // Constructor with parameters for all fields
    SVCall(uint32_t start, uint32_t end, std::string sv_type, std::string alt_allele, std::string data_type, std::string genotype, double hmm_likelihood, int support, int total_support) :
        start(start), end(end), sv_type(sv_type), alt_allele(alt_allele), data_type(data_type), genotype(genotype), hmm_likelihood(hmm_likelihood), support(support), total_support(support) {}
};

void addSVCall(std::vector<SVCall>& sv_calls, uint32_t start, uint32_t end, std::string sv_type, std::string alt_allele, std::string data_type, std::string genotype, double hmm_likelihood);

void mergeSVs(std::vector<SVCall>& sv_calls, std::unordered_map<uint32_t, uint32_t>& breakpoint_support);

void filterSVsWithLowSupport(std::vector<SVCall> &sv_calls, std::unordered_map<uint32_t, uint32_t> &breakpoint_support, int min_support);

uint32_t getSVCount(const std::vector<SVCall>& sv_calls);

void concatenateSVCalls(std::vector<SVCall>& sv_calls, const std::vector<SVCall>& sv_calls_update);

#endif // SV_OBJECT_H
