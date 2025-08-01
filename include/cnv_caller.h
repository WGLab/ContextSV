// CNVCaller: Detect CNVs and return the state sequence by SNP position
// (key = [chromosome, SNP position], value = state)

#ifndef CNV_CALLER_H
#define CNV_CALLER_H

#include "khmm.h"
#include "input_data.h"
#include "sv_types.h"
#include "sv_object.h"
#include "utils.h"

/// @cond
#include <string>
#include <vector>
#include <unordered_map>
#include <set>
// #include <mutex>
#include <shared_mutex>
#include <future>

/// @endcond

using namespace sv_types;

// SNP data is a struct containing vectors used in predicting copy number
// states. It is sorted by SNP position.
struct SNPData {
    std::vector<uint32_t> pos;
    std::vector<double> pfb;
    std::vector<double> baf;
    std::vector<double> log2_cov;
    std::vector<int> state_sequence;  // Predicted copy number states
    std::vector<bool> is_snp;  // True if the position is a SNP
    double mean_chr_cov = 0;  // Mean coverage for the chromosome

    SNPData():
        pos({}),\
        pfb({}), \
        baf({}), \
        log2_cov({}), \
        state_sequence({}), \
        is_snp({}), \
        mean_chr_cov(0) {}
};

// CNVCaller: Detect CNVs and return the state sequence by SNP position
class CNVCaller {
    private:
        std::shared_mutex& shared_mutex;

        void updateSNPData(SNPData& snp_data, uint32_t pos, double pfb, double baf, double log2_cov, bool is_snp);

        void runViterbi(const CHMM& hmm, SNPData& snp_data, std::pair<std::vector<int>, double>& prediction) const;

        // Query a region for SNPs and return the SNP data
        void querySNPRegion(std::string chr, uint32_t start_pos, uint32_t end_pos, const std::vector<uint32_t>& pos_depth_map, double mean_chr_cov, SNPData& snp_data, const InputData& input_data) const;

        // Split a region into chunks for parallel processing
        std::vector<std::string> splitRegionIntoChunks(std::string chr, uint32_t start_pos, uint32_t end_pos, int chunk_count) const;

    public:
	    CNVCaller(std::shared_mutex& shared_mutex) : shared_mutex(shared_mutex) {}

        // Define a map of CNV genotypes by HMM predicted state.
        // We only use the first 3 genotypes (0/0, 0/1, 1/1) for the VCF output.
        // Each of the 6 state predictions corresponds to a copy number state
        // (0=No predicted state)
        // 0: Unknown (No predicted state)
        // 1: 1/1 (Two copy loss: homozygous deletion, GT: 1/1 for homozygous variant)
        // 2: 0/1 (One copy loss: heterozygous deletion, GT: 0/1)
        // 3: 0/0 (Normal diploid: no copy number change, GT: 0/0 for homozygous reference)
        // 4: 1/1 (Copy neutral LOH: no copy number change, GT: 1/1 for homozygous variant)
        // 5: 2/1 (One copy gain: heterozygous duplication, GT: 1/2->0/1)
        // 6: 2/2 (Two copy gain: homozygous duplication, GT: 2/2->1/1)
        const std::unordered_map<int, Genotype> StateGenotypeMap = {
            {0, Genotype::UNKNOWN},
            {1, Genotype::HOMOZYGOUS_ALT},
            {2, Genotype::HETEROZYGOUS},
            {3, Genotype::HOMOZYGOUS_REF},
            {4, Genotype::HOMOZYGOUS_ALT},
            {5, Genotype::HETEROZYGOUS},
            {6, Genotype::HOMOZYGOUS_ALT}
        };

        // Function to get the genotype string from the state
        inline Genotype getGenotypeFromCNState(int cn_state) const {
            // return StateGenotypeMap.at(cn_state);
            try {
                return StateGenotypeMap.at(cn_state);
            } catch (const std::out_of_range& e) {
                printError("ERROR: Invalid CN state: " + std::to_string(cn_state));
                return Genotype::UNKNOWN;
            }
        }

        // Run copy number prediction for a single SV candidate, returning the
        // likelihood, predicted CNV type, genotype, and whether SNPs were found
        std::tuple<double, SVType, Genotype, int> runCopyNumberPrediction(std::string chr, const CHMM& hmm, uint32_t start_pos, uint32_t end_pos, double mean_chr_cov, const std::vector<uint32_t>& pos_depth_map, const InputData& input_data) const;

        // Run copy number prediction for SVs meeting the minimum length threshold obtained from CIGAR strings
        void runCIGARCopyNumberPrediction(std::string chr, std::vector<SVCall>& sv_candidates, const CHMM& hmm, double mean_chr_cov, const std::vector<uint32_t>& pos_depth_map, const InputData& input_data) const;

        void calculateMeanChromosomeCoverage(const std::vector<std::string>& chromosomes, std::unordered_map<std::string, std::vector<uint32_t>>& chr_pos_depth_map, std::unordered_map<std::string, double>& chr_mean_cov_map, const std::string& bam_filepath, int thread_count) const;

        void readSNPAlleleFrequencies(std::string chr, uint32_t start_pos, uint32_t end_pos, std::vector<uint32_t>& snp_pos, std::unordered_map<uint32_t, double>& snp_baf, std::unordered_map<uint32_t, double>& snp_pfb, const InputData& input_data) const;

        // Save a TSV with B-allele frequencies, log2 ratios, and copy number predictions
        void saveSVCopyNumberToTSV(SNPData& snp_data, std::string filepath, std::string chr, uint32_t start, uint32_t end, std::string sv_type, double likelihood) const;

        void saveSVCopyNumberToJSON(SNPData& before_sv, SNPData& after_sv, SNPData& snp_data, std::string chr, uint32_t start, uint32_t end, std::string sv_type, double likelihood, const std::string& filepath) const;
};

#endif // CNV_CALLER_H
