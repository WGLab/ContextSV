//
// cnv_caller.h:
// Detect CNVs
//

#ifndef CNV_CALLER_H
#define CNV_CALLER_H

#include "khmm.h"
#include "snv_caller.h"

#include <string>
#include <vector>

struct RegionCoverage {
    bool valid  = false;
    int length  = 0;
    double mean = 0;
};

class CNVCaller {
    private:
        int window_size = 10000;  // Window size (bases) for calculating the Log R Ratio
        int align_start = -1;
        int align_end   = -1;

    public:
        CNVCaller();

        /// Detect CNVs
		std::vector<double> run(std::string input_filepath, SNVCaller snv_obj);

        /// Calculate Log R Ratios
		std::vector<double> calculateLogRRatios(std::string input_filepath);

        /// Calculate region mean coverage
        RegionCoverage getMeanCoverage(std::string input_filepath, char* chr, int start_pos=1, int end_pos=1, bool entire_chr=false);

        /// Get alignment start and stop positions
        int getAlignmentEndpoints(std::string input_filepath);
};

#endif // CNV_CALLER_H
