// CNVCaller: Detect CNVs and return the state sequence by SNP position
// (key = [chromosome, SNP position], value = state)

#ifndef CNV_CALLER_H
#define CNV_CALLER_H

#include "khmm.h"
#include "input_data.h"
#include "cnv_data.h"

/// @cond
#include <string>
#include <vector>
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

class CNVCaller {
    private:
        InputData* input_data;
        double chr_mean_coverage = 0;

    public:
        CNVCaller(InputData& input_data);

        // Detect CNVs and return the state sequence by SNP position
        // (key = [chromosome, SNP position], value = state)
		CNVData run();

        // Calculate coverage log2 ratios at SNP positions
		std::vector<double> calculateLog2RatioAtSNPS(std::vector<int64_t> snp_positions);

        // Calculate the mean chromosome coverage
        double calculateMeanChromosomeCoverage();

        // Calculate region mean coverage
        double calculateWindowLogRRatio(double mean_chr_cov, int start_pos, int end_pos);

        // Read SNP positions and BAF values from the VCF file
        SNPData readSNPBAFs(std::string filepath);

        // Read SNP population frequencies from the PFB file and return a vector
        // of population frequencies for each SNP location
        std::vector<double> getSNPPopulationFrequencies(std::vector<int64_t> snp_locations);

        // Save a BED file with predicted copy number states
        void saveToBED(SNPData& snp_data, std::string filepath);

        // Save a TSV with B-allele frequencies, log 2 ratios, and copy number predictions
        void saveToTSV(SNPData& snp_data, std::string filepath);
};

#endif // CNV_CALLER_H
