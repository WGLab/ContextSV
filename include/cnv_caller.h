// CNVCaller: Detect CNVs and return the state sequence by SNP position
// (key = [chromosome, SNP position], value = state)

#ifndef CNV_CALLER_H
#define CNV_CALLER_H

#include "khmm.h"
#include "input_data.h"
#include "cnv_data.h"
#include "sv_data.h"
#include "sv_types.h"

/// @cond
#include <string>
#include <vector>
#include <unordered_map>
#include <set>
#include <mutex>

#include "snp_info.h"
/// @endcond

using namespace sv_types;

// SNP data is a struct containing vectors used in predicting copy number
// states. It is sorted by SNP position.
struct SNPData {
    std::vector<int64_t> pos;
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

// SNPValues is a struct for storing chromosome SNP information
struct SNPValues {
    int64_t pos;  // SNP location
    double pfb;  // Population frequency of the B allele
    double log2_ratio;
    double baf;  // B-allele frequency
};

// Define the comparison operator for SNPValues (by location)
struct SNPValuesCompare {
    bool operator()(const SNPValues& a, const SNPValues& b) const {
        return a.pos < b.pos;
    }
};

// Define the data structure for SNP values
using BST = std::set<SNPValues, SNPValuesCompare>;  // Binary search tree

// Define the map of chromosome to SNP BST
using ChrToSNPMap = std::unordered_map<std::string, BST>;

// Map of chromosome to SNP data
using SNPDataMap = std::unordered_map<std::string, SNPData>;

// Map of chromosome to SNP population frequency data (position -> pfb)
using PFBMap = std::unordered_map<std::string, std::map<int, double>>;

// CNVCaller: Detect CNVs and return the state sequence by SNP position
class CNVCaller {
    private:
        InputData* input_data;
        double chr_mean_coverage = 0;

        // Map of chromosome to mean coverage
        std::unordered_map<std::string, double> chr_mean_cov;

        // Define a map of CNV genotypes by HMM predicted state.
        // We only use the first 3 genotypes (0/0, 0/1, 1/1) for the VCF output.
        // Each of the 6 state predictions corresponds to a copy number state:
        // 1: 0/0 (Two copy loss: homozygous deletion, GT: 0/0)
        // 2: 1/0 (One copy loss: heterozygous deletion, GT: 0/1)
        // 3: 1/1 (Normal diploid: no copy number change, GT: 1/1)
        // 4: 1/1 (Copy neutral LOH: no copy number change, GT: 1/1)
        // 5: 2/1 (One copy gain: heterozygous duplication, GT: 1/2->0/1)
        // 6: 2/2 (Two copy gain: homozygous duplication, GT: 2/2->1/1)
        std ::map<int, std::string> cnv_genotype_map = {
            {1, "0/0"},
            {2, "0/1"},
            {3, "1/1"},
            {4, "1/1"},
            {5, "0/1"},
            {6, "1/1"}
        };

        // Define a map of CNV types by HMM predicted state.
        std ::map<int, int> cnv_type_map = {
            {1, sv_types::DEL},
            {2, sv_types::DEL},
            {3, sv_types::UNKNOWN},
            {4, sv_types::UNKNOWN},
            {5, sv_types::DUP},
            {6, sv_types::DUP}
        };

        void updateSNPData(SNPData& snp_data, int64_t pos, double pfb, double baf, double log2_cov, bool is_snp);

        void updateSNPVectors(SNPData& snp_data, std::vector<int64_t>& pos, std::vector<double>& pfb, std::vector<double>& baf, std::vector<double>& log2_cov, std::vector<int>& state_sequence, std::vector<bool>& is_snp);

        SNPData querySNPRegion(std::string chr, int64_t start_pos, int64_t end_pos, SNPInfo& snp_info, std::unordered_map<uint64_t, int> pos_depth_map, double mean_chr_cov);

        void runCopyNumberPrediction(std::string chr, SVData& sv_calls, SNPInfo& snp_info, SNPData& snp_data, CHMM hmm, int window_size);

    public:
        CNVCaller(InputData& input_data);

        // Detect CNVs and return the state sequence by SNP position
        // (key = [chromosome, SNP position], value = state)
		void run(SVData& sv_calls);

        // Calculate the mean chromosome coverage
        double calculateMeanChromosomeCoverage(std::string chr);

        // Calculate read depths for a region
        void calculateDepthsForSNPRegion(std::string chr, int start_pos, int end_pos, std::unordered_map<uint64_t, int>& pos_depth_map);

        // Calculate the log2 ratio for a region given the read depths and mean
        // chromosome coverage
        double calculateLog2Ratio(int start_pos, int end_pos, std::unordered_map<uint64_t, int>& pos_depth_map, double mean_chr_cov);

        // Read SNP positions and BAF values from the VCF file of SNP calls
        void readSNPAlleleFrequencies(std::string filepath, SNPInfo& snp_info, bool whole_genome);

        // Read SNP positions and population frequencies from the VCF file for a
        // single chromosome
        void readSNPPopAlleleFrequencies(std::string filepath, SNPDataMap& snp_data_map);

        // Read SNP population frequencies from the PFB file and return a vector
        // of population frequencies for each SNP location
        void getSNPPopulationFrequencies(SNPInfo& snp_info);

        // Save a BED file with predicted copy number states
        void saveToBED(SNPDataMap& snp_data_map, std::string filepath);

        // Save a TSV with B-allele frequencies, log 2 ratios, and copy number predictions
        void saveToTSV(SNPData& snp_data, std::string filepath);
};

#endif // CNV_CALLER_H
