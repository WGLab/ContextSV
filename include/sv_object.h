#ifndef SV_OBJECT_H
#define SV_OBJECT_H

#include <vector>
#include <memory>
#include <string>
#include <set>
#include <stdexcept>

// Struct to represent a structural variant call
struct SVCall {
    uint32_t start;
    uint32_t end;
    std::string sv_type = "NA";
    std::string alt_allele = ".";
    std::string data_type = "NA";
    std::string genotype = "./.";
    double hmm_likelihood = 0.0;
    int support = 0;

    // Comparison operator for std::set
    bool operator<(const SVCall& other) const;

    // Constructor with parameters for all fields
    SVCall(uint32_t start, uint32_t end, std::string sv_type, std::string alt_allele, std::string data_type, std::string genotype, double hmm_likelihood, int support) :
        start(start), end(end), sv_type(sv_type), alt_allele(alt_allele), data_type(data_type), genotype(genotype), hmm_likelihood(hmm_likelihood), support(support) {}
};

void addSVCall(std::set<SVCall>& sv_calls, uint32_t start, uint32_t end, std::string sv_type, std::string alt_allele, std::string data_type, std::string genotype, double hmm_likelihood);

std::vector<std::set<SVCall>> splitSVsIntoChunks(std::set<SVCall>& sv_calls, int chunk_count);

uint32_t getSVCount(const std::set<SVCall>& sv_calls);

void concatenateSVCalls(std::set<SVCall>& sv_calls, const std::set<SVCall>& sv_calls_update);

#endif // SV_OBJECT_H
