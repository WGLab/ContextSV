
#include "cnv_caller.h"
#include "utils.h"

#include <htslib/sam.h>

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>       /* log2 */
#include <algorithm>
#include <limits>
#include <tuple>
#include <iomanip>  // Progress bar


#define BUFFER_SIZE 1024

CNVCaller::CNVCaller(InputData& input_data)
{
    this->input_data = &input_data;
}

CNVData CNVCaller::run()
{
    // Read SNP positions and BAF values from the VCF file
    std::cout << "Reading SNP positions and BAF values from the VCF file..." << std::endl;
    std::pair<std::vector<int>, std::vector<double>> bafs_by_pos = readSNPBAFs();
    std::cout << "Complete." << std::endl;

    // Extract the SNP positions and BAF values
    std::vector<int> snp_positions = bafs_by_pos.first;
    std::vector<double> baf = bafs_by_pos.second;

    // Calculate LRRs
    std::cout << "Calculating LRRs at SNP positions..." << std::endl;
    std::vector<double> lrr = calculateLogRRatiosAtSNPS(snp_positions);
    std::cout << "LRRs calculated." << std::endl;

    // Read the HMM from file
    //std::string hmm_filepath = "data/wgs.hmm";
    std::string hmm_filepath = "data/hh550.hmm";
    std::cout << "Reading HMM from file: " << hmm_filepath << std::endl;
    CHMM hmm = ReadCHMM(hmm_filepath.c_str());
    std::cout << "HMM read from file." << std::endl;

    // Set up the input variables
    int num_probes = lrr.size();
    double *lrr_ptr = lrr.data();
    double *baf_ptr = baf.data();

    // Create a double array for pop. frequency and snp distance (not used), and log probabilities
    int *snpdist = NULL;
    double *pfb = NULL;
    double *logprob = NULL;

    // Run the Viterbi algorithm
    // Each of the 6 states corresponds to a CNV call:
    // 1: 0/0 (Two copy loss)
    // 2: 1/0 (One copy loss)
    // 3: 1/1 (Normal)
    // 4: 1/1 (Copy neutral LOH)
    // 5: 2/1 (One copy gain)
    // 6: 2/2 (Two copy gain)
    std::cout << "Running the Viterbi algorithm..." << std::endl;
    std::vector<int> state_sequence;  // Create the output state sequence
    state_sequence = testVit_CHMM(hmm, num_probes, lrr_ptr, baf_ptr, pfb, snpdist, logprob);
    std::cout << "Viterbi algorithm complete." << std::endl;

    // Save a CSV of the positions, LRRs, BAFs, and state sequence
    std::cout << "Saving CSV of positions, LRRs, and BAFs..." << std::endl;
    std::string output_filepath = this->input_data->getOutputDir() + "/snp_lrr_baf.csv";
    saveSNPLRRBAFCSV(output_filepath, snp_positions, baf, lrr, state_sequence);
    std::cout << "CSV saved to: " << output_filepath << std::endl;

    // Return a map of the state sequence by position
    std::string chr = this->input_data->getRegionChr();

    // Convert constant char * to char *
    char *chr_cstr = new char[chr.length() + 1];
    strcpy(chr_cstr, chr.c_str());

    // Create a map of the state sequence by position
    CNVData state_sequence_by_pos;
    //std::map<std::pair<char *, int>, int> state_sequence_by_pos;
    for (int i = 0; i < num_probes; i++)
    {
        state_sequence_by_pos.addCNVCall(chr_cstr, snp_positions[i], state_sequence[i]);
        //std::pair<char *, int> pos = std::make_pair(chr, snp_positions[i]);
        //state_sequence_by_pos[pos] = state_sequence[i];
    }

    // Free the memory
    delete[] chr_cstr;

    return state_sequence_by_pos;
}

std::vector<double> CNVCaller::calculateLogRRatiosAtSNPS(std::vector<int> snp_positions)
{
    // Get the target chromosome
    std::string chr = this->input_data->getRegionChr();
    
    // Calculate mean chromosome coverage
    std::string input_filepath = this->input_data->getBAMFilepath();
    std::cout <<  "\nCalculating coverage for chromosome: " << chr << std::endl;
    //double mean_chr_cov = calculateMeanChromosomeCoverage();  // Commented out for testing
    double mean_chr_cov = 39.4096;  // Chr6 mean coverage from test data
    // double mean_chr_cov = 39.561;  // Chr3 mean coverage from test data

    std::cout << "Mean coverage: " << mean_chr_cov << std::endl;

    // Set the region start and end from the first and last SNPs
    int region_start = snp_positions.front();
    int region_end = snp_positions.back();
    if (this->input_data->getRegionSet()) {
        region_start = std::max(region_start, this->input_data->getRegionStart());
        region_end = std::min(region_end, this->input_data->getRegionEnd());
    }

    std::cout << "Beginning analysis of region: " << chr << ":" << region_start << "-" << region_end << std::endl;

    // Loop through each SNP and calculate the LRR
    std::vector<double> snp_lrr;
    int window_size = this->input_data->getWindowSize();
    int snp_count = (int) snp_positions.size();
    for (int i = 0; i < snp_count; i++) {
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

        // Update the progress bar
        printProgress(i, snp_count-1);
        //this->input_data->printProgress(i, snp_positions.size()-1);
    }

    return snp_lrr;
}

/// Calculate the mean chromosome coverage
double CNVCaller::calculateMeanChromosomeCoverage()
{
    std::string chr = this->input_data->getRegionChr();
    std::string input_filepath = this->input_data->getBAMFilepath();

    char cmd[BUFFER_SIZE];
    FILE *fp;
    char line[BUFFER_SIZE];
    RegionCoverage cov;

    // Open a SAMtools process to calculate cumulative read depth and position
    // counts (non-zero depths only) for a single chromosome

    // Run the entire chromosome
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
    double mean_chr_cov = -1;
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
    std::string chr = this->input_data->getRegionChr();
    std::string input_filepath = this->input_data->getBAMFilepath();

    char cmd[BUFFER_SIZE];
    FILE *fp;
    char line[BUFFER_SIZE];

    // Open a SAMtools process to calculate cumulative read depth and position
    // counts (non-zero depths only) for a single region

    // Run the entire chromosome
    snprintf(cmd, BUFFER_SIZE,\
    "samtools depth -r %s:%d-%d %s | awk '{c++;s+=$3}END{print c, s}'",\
        chr.c_str(), start_pos, end_pos, input_filepath.c_str());
    //fprintf(stdout, "%s\n", cmd);  // Print the command
    //fflush(stdout);
    
    fp = popen(cmd, "r");
    if (fp == NULL)
    {
        fprintf(stderr, "Failed to run command\n");
        exit(EXIT_FAILURE);
    }

    // Parse the outputs
    uint64_t pos_count, cum_depth;
    double region_lrr = -1;
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
    std::string snp_filepath = this->input_data->getSNPFilepath();

    // Create a VCF filepath of filtered SNPs
    std::string filtered_snp_vcf_filepath = this->input_data->getOutputDir() + "/filtered_snps.vcf";

    std::cout << "Parsing SNPs from " << snp_filepath << std::endl;

    // Filter variants by depth and quality and SNPs only
    std::string region = this->input_data->getRegion();
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

void CNVCaller::saveSNPLRRBAFCSV(std::string filepath, std::vector<int> snp_positions, std::vector<double> bafs, std::vector<double> logr_ratios, std::vector<int> state_sequence)
{
    // Open the CSV file for writing
    std::ofstream csv_file(filepath);

    // Write the header
    csv_file << "position,baf,log2_ratio,cnv_state" << std::endl;

    // Write the data
    int snp_count = (int) snp_positions.size();
    for (int i = 0; i < snp_count; i++)
    {
        csv_file << snp_positions[i] << "," << bafs[i] << "," << logr_ratios[i] << "," << state_sequence[i] << std::endl;
    }

    // Close the file
    csv_file.close();
}
