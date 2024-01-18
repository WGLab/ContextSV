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
using PFBMap = std::unordered_map<std::string, std::map<int, double>>;

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
        void setAlleleFreqFilepaths(std::string filepath);
        std::string getAlleleFreqFilepaths();

        // Get the population frequency map.
        PFBMap getPFBMap();

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
        void setMeanChromosomeCoverage(std::string chr_cov);
        double getMeanChromosomeCoverage(std::string chr);

        // Set the output directory where the results will be written.
        void setOutputDir(std::string dirpath);
        std::string getOutputDir();

        // Set the number of threads to use when parallelization is possible.
        void setThreadCount(int thread_count);
        int getThreadCount();

        // Disable CIGAR string SV calling. This is useful for testing.
        void setDisableCIGAR(bool disable_cigar);
        bool getDisableCIGAR();

        // Disable SNP-based CNV calling. This is useful for calling SVs on
        // assemblies since alignment calling is sufficient.
        void setDisableSNPCNV(bool disable_snp_cnv);
        bool getDisableSNPCNV();

        // Set the filepath to the TSV file with the CNV predictions.
        void setCNVFilepath(std::string filepath);
        std::string getCNVFilepath();

        // Set the whole genome flag to true if the entire genome is being
        // analyzed.
        void setWholeGenome(bool whole_genome);
        bool getWholeGenome();

        // Read a VCF file and store the population frequencies in a map.
        void readChromosomeAFs(std::string chr, std::string filepath, std::mutex &pfb_mtx, std::mutex &print_mtx);

        // Add a chromosome's population frequency data to the map in a thread-safe manner.
        void addChromosomePopulationFrequency(std::string chr, std::map<int, double> pfb_map, std::mutex &mutex);

        // Set the verbose flag to true if verbose output is desired.
        void setVerbose(bool verbose);
        bool getVerbose();
        
    private:
        std::string short_read_bam;
        std::string long_read_bam;
        std::string ref_filepath;
        std::string snp_vcf_filepath;
        std::string pfb_filepath;
        PFBMap pfb_map;  // Map of population frequencies by SNP position (chr -> pos -> pfb)
        FASTAQuery fasta_query;
        std::string output_dir;
        std::string region;
        int window_size;
        std::string region_chr;
        int64_t region_start;
        int64_t region_end;
        bool region_set;
        std::map<std::string, double> chr_cov;  // Map of pre-calculated mean coverage values for each chromosome
        int thread_count;
        std::string hmm_filepath;
        bool disable_cigar;
        bool disable_snp_cnv;
        std::string cnv_filepath;
        bool whole_genome;  // True if the entire genome is being analyzed
        bool verbose;  // True if verbose output is enabled
};

#endif // INPUT_DATA_H
