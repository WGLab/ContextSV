//
// common.h:
// Manage common parameters

#ifndef CONTEXTSV_COMMON_H
#define CONTEXTSV_COMMON_H

#include <string>
#include <vector>
class Common {
    public:
        std::string get_bam_filepath();
        void set_bam_filepath(std::string filepath);
        std::string get_ref_filepath();
		void set_ref_filepath(std::string filepath);
        std::string get_output_dir();
		void set_output_dir(std::string dirpath);
        std::string get_region();
		void set_region(std::string region);
        int get_window_size();
		void set_window_size(int window_size);
        std::string get_snp_vcf_filepath();
        void set_snp_vcf_filepath(std::string filepath);
        std::string get_region_chr();
        int get_region_start();
        int get_region_end();

    private:
        std::string bam_filepath = "";
        std::string ref_filepath = "";
        std::string output_dir   = "";
        std::string region = "";
        int window_size = 10000;
        std::string snp_vcf_filepath = "";
        std::string region_chr = "";
        int region_start = -1;
        int region_end   = -1;
};

#endif //CONTEXTSV_COMMON_H