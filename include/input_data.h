//
// common.h:
// Manage common types, parameters, and functions

#ifndef INPUT_DATA_H
#define INPUT_DATA_H

#include <string>
#include <vector>
#include <map>

class InputData {
    public:
        std::string getBAMFilepath();
        void setBAMFilepath(std::string filepath);
        std::string getRefFilepath();
		void setRefFilepath(std::string filepath);
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

    private:
        std::string bam_filepath = "";
        std::string ref_filepath = "";
        std::string output_dir   = "";
        std::string region = "";
        int window_size = 10000;
        std::string snp_vcf_filepath = "";
        std::string region_chr = "";
        int region_start = 0;
        int region_end   = 0;
        bool region_set = false;
};

#endif // INPUT_DATA_H
