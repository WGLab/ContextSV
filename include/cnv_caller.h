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
    std::vector<std::pair<int, double>> baf_by_pos;
};

class CNVCaller {
    private:
        std::string bam_filepath;
        std::string ref_filepath;
        std::string output_dir;
        std::string region;
        std::string region_chr;
        int region_start = -1;
        int region_end   = -1;
        int window_size = 10000;  // Window size (bases) for calculating the Log R Ratio
        int align_start = -1;
        int align_end   = -1;

    public:
        CNVCaller();

        /// Detect CNVs
		std::vector<double> run(std::string region_chr, int region_start, int region_end, int window_size);

        /// Calculate Log R Ratios
		std::vector<double> calculateLogRRatios(std::string input_filepath);

        /// Calculate the mean chromosome coverage
        RegionCoverage getChromosomeCoverage(std::string input_filepath, std::string chr);

        /// Calculate region mean coverage
        RegionCoverage getRegionCoverage(std::string input_filepath, std::string chr, int start_pos=1, int end_pos=1);

        /// Get alignment start and stop positions
        int getAlignmentEndpoints(std::string input_filepath);

        // Set the bam file path
        void set_bam_filepath(std::string bam_filepath);

        // Set the reference file path
        void set_ref_filepath(std::string ref_filepath);

        // Set the output directory
        void set_output_dir(std::string output_dir);

        // Set the region to analyze
        void set_region(std::string region);
};

#endif // CNV_CALLER_H
