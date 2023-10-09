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
        std::string getPFBFilepath();
        void setPFBFilepath(std::string filepath);
        void setThreadCount(int thread_count);
        int getThreadCount();
        std::string getHMMFilepath();
        void setHMMFilepath(std::string filepath);

    private:
        std::string bam_filepath;
        std::string ref_filepath;
        std::string snp_vcf_filepath;
        std::string pfb_filepath;
        FASTAQuery fasta_query;
        std::string output_dir;
        std::string region;
        int window_size;
        std::string region_chr;
        int region_start;
        int region_end;
        bool region_set;
        std::map<std::string, double> chr_cov;
        int thread_count;
        std::string hmm_filepath;
};

#endif // INPUT_DATA_H
