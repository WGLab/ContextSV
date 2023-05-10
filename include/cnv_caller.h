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
    double baf  = 0;
};

class CNVCaller {
    private:
        int window_size = 10000;  // Window size (bases) for calculating the Log R Ratio
        int align_start = -1;
        int align_end   = -1;
        bool uses_chr_prefix = false;  // Does the BAM file use chr prefix notation?

    public:
        CNVCaller();

        /// Detect CNVs
		std::vector<double> run(std::string input_filepath, SNVCaller snv_obj);

        /// Calculate Log R Ratios
		std::vector<double> calculateLogRRatios(std::string input_filepath);

        /// Calculate the mean chromosome coverage
        RegionCoverage getChromosomeCoverage(std::string input_filepath, char *chr);

        /// Calculate region mean coverage
        RegionCoverage getRegionCoverage(std::string input_filepath, char* chr, int start_pos=1, int end_pos=1);

        /// Get alignment start and stop positions
        int getAlignmentEndpoints(std::string input_filepath);

        /// Set the chromosome prefix notation
        void setChrPrefix(bool uses_chr_prefix);
};

#endif // CNV_CALLER_H
