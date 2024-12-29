// CNVCaller: Detect CNVs and return the state sequence by SNP position
// (key = [chromosome, SNP position], value = state)

#ifndef CNV_CALLER_H
#define CNV_CALLER_H

#include "khmm.h"
#include "input_data.h"
#include "sv_types.h"
#include "sv_object.h"

/// @cond
#include <string>
#include <vector>
#include <unordered_map>
#include <set>
#include <mutex>
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
        //mutable std::mutex snp_file_mtx;  // SNP file mutex
        //mutable std::mutex pfb_file_mtx;  // Population frequency file mutex
        //mutable std::mutex bam_file_mtx;  // BAM file mutex
        std::mutex& shared_mutex;

        // Define a map of CNV genotypes by HMM predicted state.
        // We only use the first 3 genotypes (0/0, 0/1, 1/1) for the VCF output.
        // Each of the 6 state predictions corresponds to a copy number state
        // (0=No predicted state)
        // 0: Unknown (No predicted state)
        // 1: 0/0 (Two copy loss: homozygous deletion, GT: 0/0)
        // 2: 1/0 (One copy loss: heterozygous deletion, GT: 0/1)
        // 3: 1/1 (Normal diploid: no copy number change, GT: 1/1)
        // 4: 1/1 (Copy neutral LOH: no copy number change, GT: 1/1)
        // 5: 2/1 (One copy gain: heterozygous duplication, GT: 1/2->0/1)
        // 6: 2/2 (Two copy gain: homozygous duplication, GT: 2/2->1/1)
        std ::map<int, std::string> cnv_genotype_map = {
            {0, "./."},
            {1, "0/0"},
            {2, "0/1"},
            {3, "1/1"},
            {4, "1/1"},
            {5, "0/1"},
            {6, "1/1"}
        };

        void updateSNPData(SNPData& snp_data, uint32_t pos, double pfb, double baf, double log2_cov, bool is_snp);

        void runViterbi(const CHMM& hmm, SNPData& snp_data, std::pair<std::vector<int>, double>& prediction) const;

        // Query a region for SNPs and return the SNP data
        void querySNPRegion(std::string chr, uint32_t start_pos, uint32_t end_pos, const std::vector<uint32_t>& pos_depth_map, double mean_chr_cov, SNPData& snp_data, const InputData& input_data) const;

        // Split a region into chunks for parallel processing
        std::vector<std::string> splitRegionIntoChunks(std::string chr, uint32_t start_pos, uint32_t end_pos, int chunk_count) const;

    public:
        // explicit CNVCaller(const InputData& input_data);
        // Constructor with no arguments
        //CNVCaller() = default;
	    CNVCaller(std::mutex& mtx) : shared_mutex(mtx) {}

        // Run copy number prediction for a single SV candidate, returning the
        // likelihood, predicted CNV type, genotype, and whether SNPs were found
        std::tuple<double, SVType, std::string, bool> runCopyNumberPrediction(std::string chr, const CHMM& hmm, uint32_t start_pos, uint32_t end_pos, double mean_chr_cov, const std::vector<uint32_t>& pos_depth_map, const InputData& input_data) const;

        // Run copy number prediction for SVs meeting the minimum length threshold obtained from CIGAR strings
        void runCIGARCopyNumberPrediction(std::string chr, std::vector<SVCall>& sv_candidates, const CHMM& hmm, double mean_chr_cov, const std::vector<uint32_t>& pos_depth_map, const InputData& input_data) const;

        double calculateMeanChromosomeCoverage(std::string chr, std::vector<uint32_t>& chr_pos_depth_map, const std::string& bam_filepath, int thread_count, bool single_chr) const;

        void calculateRegionLog2Ratio(uint32_t start_pos, uint32_t end_pos, int sample_size, const std::vector<uint32_t>& pos_depth_map, double mean_chr_cov, std::vector<double>& pos_log2) const;

        void readSNPAlleleFrequencies(std::string chr, uint32_t start_pos, uint32_t end_pos, std::vector<uint32_t>& snp_pos, std::vector<double>& snp_baf, std::vector<double>& snp_pfb, std::vector<bool>& is_snp, const InputData& input_data) const;

        // Save a TSV with B-allele frequencies, log2 ratios, and copy number predictions
        void saveSVCopyNumberToTSV(SNPData& snp_data, std::string filepath, std::string chr, uint32_t start, uint32_t end, std::string sv_type, double likelihood) const;
};

#endif // CNV_CALLER_H
