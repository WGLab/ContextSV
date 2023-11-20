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

        std::string getShortReadBam();

        void setShortReadBam(std::string filepath);

        std::string getLongReadBam();

        void setLongReadBam(std::string filepath);

        // Set the filepath to the HMM parameters.
        void setHMMFilepath(std::string filepath);
        std::string getHMMFilepath();

        // Set the filepath to the reference genome FASTA file.
		void setRefGenome(std::string fasta_filepath);
        FASTAQuery getRefGenome();

        // Set the filepath to the tab-delimited file with SNP population frequencies.
        void setPFBFilepath(std::string filepath);
        std::string getPFBFilepath();
        
        // Set the filepath to the VCF file with SNP calls used for CNV
        // detection with the HMM.
        void setSNPFilepath(std::string filepath);
        std::string getSNPFilepath();

        // Set the genomic region to analyze.
		void setRegion(std::string region);
        std::string getRegion();
        std::string getRegionChr();
        int getRegionStart();
        int getRegionEnd();
        bool getRegionSet();

        // Set the window size for the log2 ratio calculation.
        void setWindowSize(int window_size);
        int getWindowSize();

        // Set entire-chromosome mean coverage values to speed up the log2 ratio calculations.
        void setChrCov(std::string chr_cov);
        int getChrCov(std::string chr, double& cov);

        // Set the output directory where the results will be written.
        void setOutputDir(std::string dirpath);
        std::string getOutputDir();

        // Set the number of threads to use when parallelization is possible.
        void setThreadCount(int thread_count);
        int getThreadCount();

        // Disable CIGAR string SV calling. This is useful for testing.
        void setDisableCIGAR(bool disable_cigar);
        bool getDisableCIGAR();

        // Set the filepath to the TSV file with the CNV predictions.
        void setCNVFilepath(std::string filepath);
        std::string getCNVFilepath();
        
    private:
        std::string short_read_bam;
        std::string long_read_bam;
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
        bool disable_cigar;
        std::string cnv_filepath;
};

#endif // INPUT_DATA_H
