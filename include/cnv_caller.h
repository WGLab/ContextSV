//
// cnv_caller.h:
// Detect CNVs
//

#ifndef CNV_CALLER_H
#define CNV_CALLER_H

#include "khmm.h"
#include "snv_caller.h"
#include "common.h"

#include <string>
#include <vector>

struct RegionCoverage {
    bool valid  = false;
    int length  = 0;
    double mean = 0;
    double baf  = 0;
    std::vector<std::pair<int, double>> baf_by_pos;
};

class CNVCaller {
    private:
        // std::string bam_filepath;
        // std::string ref_filepath;
        // std::string output_dir;
        // std::string region;
        // std::string region_chr;
        // int region_start = -1;
        // int region_end   = -1;
        // int window_size = 10000;  // Window size (bases) for calculating the
        // Log R Ratio
        Common common;
        int align_start = -1;
        int align_end   = -1;

    public:
        CNVCaller(Common common);

        /// Detect CNVs
		std::vector<double> run();

        /// Calculate Log R Ratios
		std::vector<double> calculateLogRRatios();

        /// Calculate the mean chromosome coverage
        RegionCoverage getChromosomeCoverage();

        /// Calculate region mean coverage
        RegionCoverage getRegionCoverage(int start_pos=1, int end_pos=1);

        /// Get alignment start and stop positions
        int getAlignmentEndpoints();
};

#endif // CNV_CALLER_H
