//
// integrative_caller.h:
// Integrate SV calls
//

#ifndef CONTEXTSV_INTEGRATIVE_CALLER_H
#define CONTEXTSV_INTEGRATIVE_CALLER_H


#include <string>
class IntegrativeCaller {
	private:
		std::string bam_filepath;
		std::string ref_filepath;
        std::string output_dir;
        std::string region;
        std::string region_chr;
        int region_start = -1;
        int region_end   = -1;
		int window_size = 10000;  // Window size (bases) for calculating the Log R Ratio

	public:
		IntegrativeCaller();

		// Entry point
		int run();

		// Check if the bam file uses chr prefix notation
		int bamHasChrPrefix(std::string filepath, bool& uses_chr_prefix);

		// Set the bam file path
		void set_bam_filepath(std::string bam_filepath);

		// Set the reference file path
		void set_ref_filepath(std::string ref_filepath);
        void set_output_dir(std::string output_dir);
        void set_region(std::string region);
};

#endif
