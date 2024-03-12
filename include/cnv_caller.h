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
#include <future>

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

// CNVCaller: Detect CNVs and return the state sequence by SNP position
class CNVCaller {
    private:
        InputData* input_data;
        double chr_mean_coverage = 0;

        // Mutex for locking the SV candidate map
        mutable std::mutex sv_candidates_mtx;

        // Mutex for locking the SNP info
        mutable std::mutex snp_data_mtx;

        // Mutex for locking the SNP data
        mutable std::mutex snp_mtx;

        // Mutex for locking the position-depth map
        mutable std::mutex pos_depth_map_mtx;

        // Mutex for locking the HMM prediction
        mutable std::mutex hmm_mtx;

        // Map of chromosome to mean coverage
        std::unordered_map<std::string, double> chr_mean_cov;

        // Define a map of CNV genotypes by HMM predicted state.
        // We only use the first 3 genotypes (0/0, 0/1, 1/1) for the VCF output.
        // Each of the 6 state predictions corresponds to a copy number state
        // (0=No predicted state)
        // 0: 1/1 (Normal diploid: no copy number change, GT: 1/1)
        // 1: 0/0 (Two copy loss: homozygous deletion, GT: 0/0)
        // 2: 1/0 (One copy loss: heterozygous deletion, GT: 0/1)
        // 3: 1/1 (Normal diploid: no copy number change, GT: 1/1)
        // 4: 1/1 (Copy neutral LOH: no copy number change, GT: 1/1)
        // 5: 2/1 (One copy gain: heterozygous duplication, GT: 1/2->0/1)
        // 6: 2/2 (Two copy gain: homozygous duplication, GT: 2/2->1/1)
        std ::map<int, std::string> cnv_genotype_map = {
            {0, "1/1"},
            {1, "0/0"},
            {2, "0/1"},
            {3, "1/1"},
            {4, "1/1"},
            {5, "0/1"},
            {6, "1/1"}
        };

        // Define a map of CNV types by HMM predicted state (0=No predicted state)
        std ::map<int, int> cnv_type_map = {
            {0, sv_types::UNKNOWN},
            {1, sv_types::DEL},
            {2, sv_types::DEL},
            {3, sv_types::UNKNOWN},
            {4, sv_types::UNKNOWN},
            {5, sv_types::DUP},
            {6, sv_types::DUP}
        };

        void updateSNPData(SNPData& snp_data, int64_t pos, double pfb, double baf, double log2_cov, bool is_snp);

        void updateSNPVectors(SNPData& snp_data, std::vector<int64_t>& pos, std::vector<double>& pfb, std::vector<double>& baf, std::vector<double>& log2_cov, std::vector<int>& state_sequence, std::vector<bool>& is_snp);

        std::vector<int> runViterbi(CHMM hmm, SNPData &snp_data);

        std::pair<SNPData, bool> querySNPRegion(std::string chr, int64_t start_pos, int64_t end_pos, SNPInfo &snp_info, std::unordered_map<uint64_t, int> &pos_depth_map, double mean_chr_cov);

        SNPData runCopyNumberPrediction(std::string chr, std::map<SVCandidate, SVInfo>& sv_candidates, SNPInfo& snp_info, CHMM hmm, int window_size, double mean_chr_cov);

        SNPData runCopyNumberPredictionChunk(std::string chr, std::map<SVCandidate, SVInfo>& sv_candidates, std::vector<SVCandidate> sv_chunk, SNPInfo& snp_info, CHMM hmm, int window_size, double mean_chr_cov, std::unordered_map<uint64_t, int>& pos_depth_map);

        void updateSVType(std::map<SVCandidate, SVInfo>& sv_candidates, SVCandidate key, int sv_type, std::string data_type);

        void updateSVGenotype(std::map<SVCandidate, SVInfo>& sv_candidates, SVCandidate key, std::string genotype);

        std::vector<std::string> splitRegionIntoChunks(std::string chr, int64_t start_pos, int64_t end_pos, int chunk_count);

        std::vector<std::vector<SVCandidate>> splitSVCandidatesIntoChunks(std::map<SVCandidate, SVInfo>& sv_candidates, int chunk_count);

        void mergePosDepthMaps(std::unordered_map<uint64_t, int>& pos_depth_map, std::unordered_map<uint64_t, int>& pos_depth_map_chunk);

    public:
        CNVCaller(InputData& input_data);

        // Detect CNVs and return the state sequence by SNP position
        // (key = [chromosome, SNP position], value = state)
		void run(SVData& sv_calls);

        // Calculate the mean chromosome coverage
        double calculateMeanChromosomeCoverage(std::string chr);

        // Calculate read depths for a region
        void calculateDepthsForSNPRegion(std::string chr, int64_t start_pos, int64_t end_pos, std::unordered_map<uint64_t, int>& pos_depth_map);

        // Calculate the log2 ratio for a region given the read depths and mean
        // chromosome coverage
        double calculateLog2Ratio(int start_pos, int end_pos, std::unordered_map<uint64_t, int>& pos_depth_map, double mean_chr_cov);

        // Read SNP positions and BAF values from the VCF file of SNP calls
        void readSNPAlleleFrequencies(std::string chr, std::string filepath, SNPInfo& snp_info);

        // Read SNP population frequencies from the PFB file and return a vector
        // of population frequencies for each SNP location
        void getSNPPopulationFrequencies(std::string chr, SNPInfo& snp_info);

        // Save a TSV with B-allele frequencies, log 2 ratios, and copy number predictions
        void saveToTSV(SNPData& snp_data, std::string filepath);
};

#endif // CNV_CALLER_H
