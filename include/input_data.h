//
// common.h:
// Manage common types, parameters, and functions

#ifndef INPUT_DATA_H
#define INPUT_DATA_H

#include "fasta_query.h"

/// @cond
#include <string>
#include <vector>
#include <map>
/// @endcond

class InputData {
    public:
        InputData();
        std::string getBAMFilepath();
        void setBAMFilepath(std::string filepath);
		void setRefGenome(std::string fasta_filepath);
        FASTAQuery getRefGenome();
        std::string getSequence(std::string chr, int pos_start, int pos_end);
        std::string getOutputDir();
        void setOutputDir(std::string dirpath);
        std::string getRegion();
		void setRegion(std::string region);
        int getWindowSize();
		void setWindowSize(int window_size);
        std::string getSNPFilepath();
        void setSNPFilepath(std::string filepath);
        std::string getRegionChr();
        int getRegionStart();
        int getRegionEnd();
        bool getRegionSet();
        void setChrCov(std::string chr_cov);
        int getChrCov(std::string chr, double& cov);

    private:
        std::string bam_filepath = "";
        std::string ref_filepath = "";
        FASTAQuery fasta_query;
        std::string output_dir   = "";
        std::string region = "";
        int window_size = 10000;
        std::string snp_vcf_filepath = "";
        std::string region_chr = "";
        int region_start = 0;
        int region_end   = 0;
        bool region_set = false;
        std::map<std::string, double> chr_cov;
};

#endif // INPUT_DATA_H
