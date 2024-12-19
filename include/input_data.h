//
// common.h:
// Manage common types, parameters, and functions

#ifndef INPUT_DATA_H
#define INPUT_DATA_H

#include "fasta_query.h"

/// @cond
#include <string>
#include <vector>
// #include <map>
#include <map>
#include <unordered_map>
#include <mutex>
/// @endcond

// Type definition for B-allele population frequency map (chr -> pos -> pfb)
// using PFBMap = std::unordered_map<std::string, std::map<int, double>>;

class InputData {
    public:
        InputData();

        std::string getShortReadBam() const;

        void setShortReadBam(std::string filepath);

        std::string getLongReadBam() const;

        void setLongReadBam(std::string filepath);

        // Set the filepath to the HMM parameters.
        void setHMMFilepath(std::string filepath);
        std::string getHMMFilepath() const;

        // Set the filepath to the reference genome FASTA file.
		void setRefGenome(std::string filepath);
        std::string getRefGenome() const;

        // Set the filepath to the text file containing the locations of the
        // VCF files with population frequencies for each chromosome.
        void setAlleleFreqFilepaths(std::string filepath);
        std::string getAlleleFreqFilepath(std::string chr) const;

        // Set the filepath to the VCF file with SNP calls used for CNV
        // detection with the HMM.
        void setSNPFilepath(std::string filepath);
        std::string getSNPFilepath() const;

        // Set the ethnicity for SNP population frequencies.
        void setEthnicity(std::string ethnicity);
        std::string getEthnicity() const;

        // Set the sample size for HMM predictions.
        void setSampleSize(int sample_size);
        int getSampleSize() const;

        // Set the minimum CNV length to use for copy number predictions.
        void setMinCNVLength(int min_cnv_length);
        uint32_t getMinCNVLength() const;

        // Set the chromosome to analyze.
        void setChromosome(std::string chr);
        std::string getChromosome() const;
        bool isSingleChr() const;

        // Set the region to analyze.
        void setRegion(std::string region);
        std::pair<int32_t, int32_t> getRegion() const;
        bool isRegionSet() const;

        // Set the output directory where the results will be written.
        void setOutputDir(std::string dirpath);
        std::string getOutputDir() const;

        // Set the number of threads to use when parallelization is possible.
        void setThreadCount(int thread_count);
        int getThreadCount() const;

        // Set the verbose flag to true if verbose output is desired.
        void setVerbose(bool verbose);
        bool getVerbose();

        // Set whether to extend the SNP CNV regions around the SV breakpoints
        // (+/- 1/2 SV length), save a TSV file, and generate HTML reports.
        void saveCNVData(bool save_cnv_data);
        bool getSaveCNVData() const;
        
    private:
        std::string short_read_bam;
        std::string long_read_bam;
        std::string ref_filepath;
        std::string snp_vcf_filepath;
        std::string ethnicity;
        std::unordered_map<std::string, std::string> pfb_filepaths;  // Map of population frequency VCF filepaths by chromosome
        std::string output_dir;
        int sample_size;
        uint32_t min_cnv_length;
        std::string chr;  // Chromosome to analyze
        std::pair<int32_t, int32_t> start_end;  // Region to analyze
        bool region_set;  // True if a region is set
        int thread_count;
        std::string hmm_filepath;
        std::string cnv_filepath;
        bool verbose;  // True if verbose output is enabled
        bool save_cnv_data;  // True if SNP CNV regions should be extended around SV breakpoints, and saved to a TSV file (Large performance hit)
        bool single_chr;
};

#endif // INPUT_DATA_H
