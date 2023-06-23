
// #include "khmm.h"
#include "cnv_caller.h"

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>       /* log2 */
#include <htslib/sam.h>
#include <algorithm>
#include <limits>
#include <tuple>
#include "common.h"

#define BUFFER_SIZE 1024

CNVCaller::CNVCaller(Common common)
{
    this->common = common;
}

std::vector<double> CNVCaller::run()
{
    // Read SNP positions and BAF values from the VCF file
    std::cout << "Reading SNP positions and BAF values from the VCF file..." << std::endl;
    std::pair<std::vector<int>, std::vector<double>> bafs_by_pos = readSNPBAFs();
    std::cout << "Complete." << std::endl;

    // Extract the SNP positions and BAF values
    std::vector<int> snp_positions = bafs_by_pos.first;
    std::vector<double> snp_bafs = bafs_by_pos.second;

    // Calculate LRRs
    std::cout << "Calculating LRRs at SNP positions..." << std::endl;
    std::vector<double> log_r_ratios = calculateLogRRatiosAtSNPS(snp_positions);
    std::cout << "LRRs calculated." << std::endl;

    // Save a CSV of the positions, LRRs, and BAFs
    std::cout << "Saving CSV of positions, LRRs, and BAFs..." << std::endl;
    std::string output_filepath = this->common.get_output_dir() + "/snp_lrr_baf.csv";
    saveSNPLRRBAFCSV(output_filepath, snp_positions, snp_bafs, log_r_ratios);
    std::cout << "CSV saved to: " << output_filepath << std::endl;

    // Calculate BAFs
    //std::vector<double> b_allele_freqs;
    //b_allele_freqs = calculateBAFs(input_filepath);

    // Read the HMM from file
    //std::string hmm_filepath = "data/wgs.hmm";
    //CHMM hmm = ReadCHMM(hmm_filepath.c_str());
    

    // Set up the input variables
    int num_probes = log_r_ratios.size();
    // double *lrr = &log_r_ratios[0];
    double *baf = NULL;
    double *pfb = NULL;
    int *snpdist = NULL;
    double logprob = 0.0;

    // Run the Viterbi algorithm
    //testVit_CHMM(hmm, num_probes, lrr, baf, pfb, snpdist, &logprob);

    // Estimate the hidden states from the LRRs
    // TODO: Follow detect_cnv.pl's example and use the Viterbi algorithm
    // https://github.com/WGLab/PennCNV/blob/b6d76b58821deea4f6fe9dc3c241215f25a7fd67/detect_cnv.pl#LL903C20-L903C20

    // #generate CNV calls
	// 			my $probe_count = scalar (@$lrr)-1;
	// 			khmm::testVit_CHMM ($hmm_model, $probe_count, $lrr, $baf, $pfb, $snpdist, \$logprob);
	// 			analyzeStateSequence ($curcnvcall, $curchr, $pfb, $name, $pos, $sample_sex);

    
    //testVit_CHMM(hmm, log_r_ratios);


    return log_r_ratios;
}

std::vector<double> CNVCaller::calculateLogRRatiosAtSNPS(std::vector<int> snp_positions)
{
    std::string target_chr = this->common.get_region_chr();
    
    // Calculate mean chromosome coverage
    std::string input_filepath = this->common.get_bam_filepath();
    std::cout <<  "\nCalculating coverage for chromosome: \n" << target_chr.c_str() << std::endl;
    double mean_chr_cov = calculateMeanChromosomeCoverage();

    std::cout << "Mean coverage: " << mean_chr_cov << std::endl;

    // Set the region start and end from the first and last SNPs
    int region_start = snp_positions.front();
    int region_end = snp_positions.back();
    if (this->common.get_region_set()) {
        region_start = std::max(region_start, this->common.get_region_start());
        region_end = std::min(region_end, this->common.get_region_end());
    }

    std::cout << "Beginning analysis of region: " << target_chr << ":" << region_start << "-" << region_end << std::endl;

    // Loop through each SNP and calculate the LRR
    std::vector<double> snp_lrr;
    int window_size = this->common.get_window_size();
    for (int i = 0; i < snp_positions.size(); i++) {
        int pos = snp_positions[i];

        // Skip SNPs outside of the region
        if (pos < region_start || pos > region_end) {
            continue;
        }

        // Calculate window mean coverage
        int window_start = pos - (window_size / 2);
        int window_end = pos + (window_size / 2);
        double lrr = calculateWindowLogRRatio(mean_chr_cov, window_start, window_end);

        // Set the LRR value
        snp_lrr.push_back(lrr);
    }

    return snp_lrr;
}

/// Calculate the mean chromosome coverage
double CNVCaller::calculateMeanChromosomeCoverage()
{
    std::string chr = this->common.get_region_chr();
    std::string input_filepath = this->common.get_bam_filepath();

    char cmd[BUFFER_SIZE];
    FILE *fp;
    char line[BUFFER_SIZE];
    RegionCoverage cov;
    bool log_debug = false;  // Log debugging output

    // Open a SAMtools process to calculate cumulative read depth and position
    // counts (non-zero depths only) for a single chromosome

    // Run the entire chromosome
    log_debug = true;
    snprintf(cmd, BUFFER_SIZE,\
    "samtools depth -r %s %s | awk '{c++;s+=$3}END{print c, s}'",\
    chr.c_str(), input_filepath.c_str());  // Remove '-a' for debugging

    // Parse the output
    fp = popen(cmd, "r");
    if (fp == NULL) {
        fprintf(stderr, "Failed to run command\n");
        exit(1);
    }
    
    fprintf(stdout, "%s\n", cmd);  // Print the command
    fflush(stdout);
    
    fp = popen(cmd, "r");
    if (fp == NULL)
    {
        fprintf(stderr, "Failed to run command\n");
        exit(EXIT_FAILURE);
    }

    // Parse the outputs
    uint64_t pos_count, cum_depth;
    double mean_chr_cov;
    if (fgets(line, BUFFER_SIZE, fp) != NULL)
    {           
        if (sscanf(line, "%ld%ld", &pos_count, &cum_depth) == 2)
        {
            // Calculate the mean chromosome coverage
            mean_chr_cov = (double) cum_depth / (double) pos_count;
        }
    }
    pclose(fp);  // Close the process

    return mean_chr_cov;
}

double CNVCaller::calculateWindowLogRRatio(double mean_chr_cov, int start_pos, int end_pos)
{
    std::string chr = this->common.get_region_chr();
    std::string input_filepath = this->common.get_bam_filepath();

    char cmd[BUFFER_SIZE];
    FILE *fp;
    char line[BUFFER_SIZE];

    // Open a SAMtools process to calculate cumulative read depth and position
    // counts (non-zero depths only) for a single region

    // Run the entire chromosome
    snprintf(cmd, BUFFER_SIZE,\
    "samtools depth -r %s:%d-%d %s | awk '{c++;s+=$3}END{print c, s}'",\
        chr.c_str(), start_pos, end_pos, input_filepath.c_str());
    fprintf(stdout, "%s\n", cmd);  // Print the command
    fflush(stdout);
    
    fp = popen(cmd, "r");
    if (fp == NULL)
    {
        fprintf(stderr, "Failed to run command\n");
        exit(EXIT_FAILURE);
    }

    // Parse the outputs
    uint64_t pos_count, cum_depth;
    double region_lrr;
    if (fgets(line, BUFFER_SIZE, fp) != NULL)
    {           
        if (sscanf(line, "%ld%ld", &pos_count, &cum_depth) == 2)
        {
            // Calculate the LRR
            double mean_window_cov = (double) cum_depth / (double) pos_count;
            region_lrr = log2(mean_window_cov / mean_chr_cov);
        }
    }
    pclose(fp);  // Close the process

    return region_lrr;
}

std::pair<std::vector<int>, std::vector<double>> CNVCaller::readSNPBAFs()
{
    // Get the VCF filepath with SNPs
    std::string snp_filepath = this->common.get_snp_vcf_filepath();

    // Create a VCF filepath of filtered SNPs
    std::string filtered_snp_vcf_filepath = this->common.get_output_dir() + "/filtered_snps.vcf";

    std::cout << "Parsing SNPs from " << snp_filepath << std::endl;

    // Filter variants by depth and quality and SNPs only
    std::string region = this->common.get_region();
    std::string cmd;
    if (region == "")
    {
        cmd = "bcftools view -v snps -i 'QUAL > 30 && DP > 10 && FILTER = \"PASS\"' " + snp_filepath + " > " + filtered_snp_vcf_filepath;
    } else {
        cmd = "bcftools view -r " + region + " -v snps -i 'QUAL > 30 && DP > 10 && FILTER = \"PASS\"' " + snp_filepath + " > " + filtered_snp_vcf_filepath;
    }
    std::cout << "Command: " << cmd << std::endl;
    system(cmd.c_str());

    std::cout << "Filtered SNPs written to " << filtered_snp_vcf_filepath << std::endl;

    // Extract all BAFs from the filtered SNPs and store in the pair vector
    std::cout << "Extracting BAFs from filtered SNPs" << std::endl;
    cmd = "bcftools query -f '%CHROM,%POS,[%VAF]\n' " + filtered_snp_vcf_filepath;
    std::cout << "Command: " << cmd << std::endl;
    FILE *fp = popen(cmd.c_str(), "r");
    if (fp == NULL)
    {
        std::cerr << "ERROR: Could not open pipe for command: " << cmd << std::endl;
        exit(1);
    }

    // Read the BAFs
    char line[BUFFER_SIZE];

    // Loop through the lines
    std::vector<int> snp_positions;
    std::vector<double> snp_bafs;
    while (fgets(line, BUFFER_SIZE, fp) != NULL)
    {
        // Parse the line
        char *tok = strtok(line, ",");  // Tokenize the line
        int col = 0;  // Column index
        std::string chr = "";
        uint64_t pos = 0;
        double baf = 0.0;
        while (tok != NULL)
        {
            // Get the chromosome from column 1
            if (col == 0)
            {
                chr = tok;
            }

            // Get the position from column 2
            else if (col == 1)
            {
                pos = atoi(tok);
            }

            // Get the BAF from column 3
            else if (col == 2)
            {
                baf = atof(tok);
            }

            tok = strtok(NULL, ",");
            col++;
        }

        // Store the position and BAF
        snp_positions.push_back(pos);
        snp_bafs.push_back(baf);
    }

    // Close the pipe
    pclose(fp);

    // Create the pair vector
    std::pair<std::vector<int>, std::vector<double>> snp_data = std::make_pair(snp_positions, snp_bafs);

    return snp_data;
}

void CNVCaller::saveSNPLRRBAFCSV(std::string filepath, std::vector<int> snp_positions, std::vector<double> bafs, std::vector<double> logr_ratios)
{
    // Open the CSV file for writing
    std::ofstream csv_file(filepath);

    // Write the header
    csv_file << "position,baf,log2_ratio" << std::endl;

    // Write the data
    for (int i = 0; i < snp_positions.size(); i++)
    {
        csv_file << snp_positions[i] << "," << bafs[i] << "," << logr_ratios[i] << std::endl;
    }

    // Close the file
    csv_file.close();
}
