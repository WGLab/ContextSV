// CNVCaller: Detect CNVs and return the state sequence by SNP position
// (key = [chromosome, SNP position], value = state)

#ifndef CNV_CALLER_H
#define CNV_CALLER_H

#include "khmm.h"
#include "input_data.h"
#include "cnv_data.h"
#include "region.h"

/// @cond
#include <string>
#include <vector>
#include <unordered_map>
/// @endcond

// SNP data is a struct containing vectors used in predicting copy number states
struct SNPData {
    std::vector<int64_t> locations;
    std::vector<double> pfbs;
    std::vector<double> bafs;
    std::vector<double> log2_ratios;
    std::vector<int> state_sequence;

    SNPData():
        locations({}),\
        pfbs({}), \
        bafs({}), \
        log2_ratios({}), \
        state_sequence({}) {}
};

// Map of chromosome to SNP data
using SNPDataMap = std::unordered_map<std::string, SNPData>;

// CNVCaller: Detect CNVs and return the state sequence by SNP position
class CNVCaller {
    private:
        InputData* input_data;
        double chr_mean_coverage = 0;

    public:
        CNVCaller(InputData& input_data);

        // Detect CNVs and return the state sequence by SNP position
        // (key = [chromosome, SNP position], value = state)
		void run(CNVData& cnv_calls);

        // Calculate coverage log2 ratios at SNP positions
		void calculateLog2RatioAtSNPS(SNPDataMap& snp_data_map);

        // Calculate the mean chromosome coverage
        double calculateMeanChromosomeCoverage(std::string chr);

        // Calculate region mean coverage
        double calculateWindowLogRRatio(double mean_chr_cov, std::string chr, int start_pos, int end_pos);

        // Read SNP positions and BAF values from the VCF file
        void readSNPAlleleFrequencies(std::string filepath, SNPDataMap& snp_data_map, bool whole_genome);

        // Read SNP population frequencies from the PFB file and return a vector
        // of population frequencies for each SNP location
        void getSNPPopulationFrequencies(PFBMap& pfb_map, SNPDataMap& snp_data_map);

        // Save a BED file with predicted copy number states
        void saveToBED(SNPDataMap& snp_data_map, std::string filepath);

        // Save a TSV with B-allele frequencies, log 2 ratios, and copy number predictions
        void saveToTSV(SNPDataMap& snp_data_map, std::string filepath);
};

#endif // CNV_CALLER_H
