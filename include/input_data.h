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

        std::string getShortReadBam();

        void setShortReadBam(std::string filepath);

        std::string getLongReadBam();

        void setLongReadBam(std::string filepath);

        // Set the filepath to the HMM parameters.
        void setHMMFilepath(std::string filepath);
        std::string getHMMFilepath();

        // Set the filepath to the reference genome FASTA file.
		void setRefGenome(std::string fasta_filepath);

        // Return a reference to the ReferenceGenome object.
        const ReferenceGenome& getRefGenome() const;
        // FASTAQuery getRefGenome();

        // Query the reference genome for a sequence.
        std::string queryRefGenome(const std::string& chr, uint32_t pos_start, uint32_t pos_end) const;

        // Get the chromosomes in the reference genome.
        std::vector<std::string> getRefGenomeChromosomes();

        // Get a chromosome's length in the reference genome.
        int64_t getRefGenomeChromosomeLength(std::string chr);

        // Set the filepath to the text file containing the locations of the
        // VCF files with population frequencies for each chromosome.
        void setAlleleFreqFilepaths(std::string filepath);

        // Get the chromosome's VCF filepath with population frequencies.
        std::string getAlleleFreqFilepath(std::string chr);

        // Get the population frequency map.
        // PFBMap getPFBMap();

        // Set the filepath to the VCF file with SNP calls used for CNV
        // detection with the HMM.
        void setSNPFilepath(std::string filepath);
        std::string getSNPFilepath();

        // Set the ethnicity for SNP population frequencies.
        void setEthnicity(std::string ethnicity);
        std::string getEthnicity();

        // Set the window size for the log2 ratio calculation.
        void setWindowSize(int window_size);
        int getWindowSize();

        // Set the minimum CNV length to use for copy number predictions.
        void setMinCNVLength(int min_cnv_length);
        int getMinCNVLength();

        // Set the chromosome to analyze.
        void setChromosome(std::string chr);
        std::string getChromosome();

        // Set the region to analyze.
        void setRegion(std::string region);
        std::pair<int32_t, int32_t> getRegion();
        bool isRegionSet();

        // Set the output directory where the results will be written.
        void setOutputDir(std::string dirpath);
        std::string getOutputDir();

        // Set the number of threads to use when parallelization is possible.
        void setThreadCount(int thread_count);
        int getThreadCount();

        // Set the verbose flag to true if verbose output is desired.
        void setVerbose(bool verbose);
        bool getVerbose();

        // Set whether to extend the SNP CNV regions around the SV breakpoints
        // (+/- 1/2 SV length), save a TSV file, and generate HTML reports.
        void saveCNVData(bool save_cnv_data);
        bool getSaveCNVData();
        
    private:
        std::string short_read_bam;
        std::string long_read_bam;
        std::string ref_filepath;
        std::string snp_vcf_filepath;
        std::string ethnicity;
        std::unordered_map<std::string, std::string> pfb_filepaths;  // Map of population frequency VCF filepaths by chromosome
        ReferenceGenome fasta_query;
        std::string output_dir;
        int window_size;
        int min_cnv_length;
        std::string chr;  // Chromosome to analyze
        std::pair<int32_t, int32_t> start_end;  // Region to analyze
        bool region_set;  // True if a region is set
        int thread_count;
        std::string hmm_filepath;
        std::string cnv_filepath;
        bool verbose;  // True if verbose output is enabled
        bool save_cnv_data;  // True if SNP CNV regions should be extended around SV breakpoints, and saved to a TSV file (Large performance hit)
};

#endif // INPUT_DATA_H
