//
// cnv_caller.h:
// Detect CNVs
//

#ifndef CNV_CALLER_H
#define CNV_CALLER_H

#include "khmm.h"
#include "common.h"

#include <string>
#include <vector>

struct RegionCoverage {
    bool valid  = false;
    uint64_t length  = 0;
    double mean = 0;
    double baf  = 0;
};

class CNVCaller {
    private:
        Common common;
        double chr_mean_coverage = 0;

    public:
        CNVCaller(Common common);

        // Detect CNVs and return the state sequence by SNP position
		std::map<int, int> run();

        // Calculate Log R Ratios
		std::vector<double> calculateLogRRatiosAtSNPS(std::vector<int> snp_positions);

        // Calculate the mean chromosome coverage
        double calculateMeanChromosomeCoverage();

        // Calculate region mean coverage
        double calculateWindowLogRRatio(double mean_chr_cov, int start_pos, int end_pos);

        // Read SNP positions and BAF values from the VCF file
        std::pair<std::vector<int>, std::vector<double>> readSNPBAFs();

        // Save a CSV with SNP positions, BAF values and Log R Ratios
        void saveSNPLRRBAFCSV(std::string filepath, std::vector<int> snp_positions, std::vector<double> bafs, std::vector<double> logr_ratios, std::vector<int> state_sequence);
};

#endif // CNV_CALLER_H
