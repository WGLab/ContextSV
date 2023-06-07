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
    uint64_t length  = 0;
    double mean = 0;
    double baf  = 0;
    std::vector<std::pair<int, double>> baf_by_pos;
};

class CNVCaller {
    private:
        Common common;
        std::vector<int> snp_positions;
        int align_start = -1;
        int align_end   = -1;

    public:
        CNVCaller(Common common, std::vector<int> snp_positions);

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
