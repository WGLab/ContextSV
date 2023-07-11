//
// integrative_caller.h:
// Integrate SV calls
//

#ifndef CONTEXTSV_INTEGRATIVE_CALLER_H
#define CONTEXTSV_INTEGRATIVE_CALLER_H

#include "common.h"

#include <string>

class IntegrativeCaller {
	private:
		Common common;

	public:
		IntegrativeCaller(Common common);

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
        void set_window_size(int window_size);
};

#endif
