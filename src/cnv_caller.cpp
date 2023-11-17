
#include "cnv_caller.h"
#include "utils.h"

#include <htslib/sam.h>

/// @cond
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
#define MIN_PFB 0.01
#define MAX_PFB 0.99
/// @endcond


CNVCaller::CNVCaller(InputData& input_data)
{
    this->input_data = &input_data;
}

CNVData CNVCaller::run()
{
    // Read SNP positions and BAF values from the VCF file
    std::cout << "Reading SNP positions and BAF values from the VCF file..." << std::endl;
    std::string snp_filepath = this->input_data->getSNPFilepath();
    SNPData snp_data = readSNPBAFs(snp_filepath);
    int snp_count = (int) snp_data.locations.size();
    std::cout << "Complete. Number of SNPs: " << snp_count << std::endl;

    // Get the population frequencies for each SNP
    if (this->input_data->getPFBFilepath() == "")
    {
        std::cout << "No PFB file provided. Using the minimum PFB value of " << MIN_PFB << " for all SNPs." << std::endl;
        snp_data.pfbs = std::vector<double>(snp_count, MIN_PFB);
    } else {
        std::cout << "Using PFB file: " << this->input_data->getPFBFilepath() << std::endl;
        std::cout << "Getting population frequencies for each SNP..." << std::endl;
        snp_data.pfbs = getSNPPopulationFrequencies(snp_data.locations);
        std::cout << "Population frequencies retrieved." << std::endl;
    }

    // Calculate LRRs
    std::cout << "Calculating log2 ratio of coverage at SNP positions..." << std::endl;
    snp_data.log2_ratios = calculateLog2RatioAtSNPS(snp_data.locations);
    std::cout << "Log2 ratios calculated." << std::endl;

    // Read the HMM from file
    //std::string hmm_filepath = "data/wgs.hmm";
    std::string hmm_filepath = this->input_data->getHMMFilepath();
    std::cout << "Reading HMM from file: " << hmm_filepath << std::endl;
    CHMM hmm = ReadCHMM(hmm_filepath.c_str());
    std::cout << "HMM read from file." << std::endl;

    // Set up the data for the Viterbi algorithm
    double *lrr_ptr = snp_data.log2_ratios.data();
    double *baf_ptr = snp_data.bafs.data();
    double *pfb_ptr = snp_data.pfbs.data();
    int *snpdist = NULL;
    double *logprob = NULL;

    // Run the Viterbi algorithm
    // Each of the 6 state predictions corresponds to a copy number state:
    // 1: 0/0 (Two copy loss)
    // 2: 1/0 (One copy loss)
    // 3: 1/1 (Normal)
    // 4: 1/1 (Copy neutral LOH)
    // 5: 2/1 (One copy gain)
    // 6: 2/2 (Two copy gain)
    std::cout << "Running the Viterbi algorithm..." << std::endl;
    std::vector<int> state_sequence;  // Create the output state sequence
    state_sequence = testVit_CHMM(hmm, snp_count, lrr_ptr, baf_ptr, pfb_ptr, snpdist, logprob);
    snp_data.state_sequence = state_sequence;
    std::cout << "Complete." << std::endl;

    // Save a TSV of the positions, LRRs, BAFs, and state sequence
    std::cout << "Saving TSV of copy number prediction data..." << std::endl;
    std::string output_tsv = this->input_data->getOutputDir() + "/cnv_data.tsv";
    saveToTSV(snp_data, output_tsv);
    std::cout << "Saved to: " << output_tsv << std::endl;

    // Save a BED file of the CNV state sequence
    std::cout << "Saving BED of CNV states..." << std::endl;
    std::string output_bed = this->input_data->getOutputDir() + "/cnv_states.bed";
    saveToBED(snp_data, output_bed);
    std::cout << "Saved to: " << output_bed << std::endl;

    // Get the chromosome as a C string
    std::string chr = this->input_data->getRegionChr();
    char *chr_cstr = new char[chr.length() + 1];
    strcpy(chr_cstr, chr.c_str());

    // Create a map of the state sequence by position
    CNVData cnv_data;
    for (int i = 0; i < snp_count; i++)
    {
        cnv_data.addCNVCall(chr_cstr, snp_data.locations[i], snp_data.state_sequence[i]);
    }

    // Free the C string
    delete[] chr_cstr;

    return cnv_data;
}

std::vector<double> CNVCaller::calculateLog2RatioAtSNPS(std::vector<int64_t> snp_locations)
{
    std::string input_filepath = this->input_data->getBAMFilepath();
    std::string chr = this->input_data->getRegionChr();

    // Check if the chromosome coverage was passed in
    double mean_chr_cov = -1;
    if (this->input_data->getChrCov(chr, mean_chr_cov) == -1)
    {
        // Calculate the mean chromosome coverage
        std::cout <<  "\nCalculating coverage for chromosome: " << chr << std::endl;
        mean_chr_cov = calculateMeanChromosomeCoverage();
        std::cout << "Mean coverage for chromosome " << chr << ": " << mean_chr_cov << std::endl;
    } else {
        std::cout << "Using user-provided mean coverage for chromosome " << chr << ": " << mean_chr_cov << std::endl;
    }

    // Set the region start and end from the first and last SNPs
    std::cout << "Getting region start and end..." << std::endl;
    int region_start = snp_locations.front();
    int region_end = snp_locations.back();
    std::cout << "Fixing region start and end..." << std::endl;
    if (this->input_data->getRegionSet()) {
        region_start = std::max(region_start, this->input_data->getRegionStart());
        region_end = std::min(region_end, this->input_data->getRegionEnd());
    }

    std::cout << "Predicting CNV states for SNPs in region: " << chr << ":" << region_start << "-" << region_end << std::endl;

    // Loop through each SNP and calculate the LRR
    std::vector<double> snp_lrr;
    int window_size = this->input_data->getWindowSize();
    int snp_count = (int) snp_locations.size();
    for (int i = 0; i < snp_count; i++) {
        int pos = snp_locations[i];

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

        // Update the progress bar every 1000 SNPs
        if (i % 1000 == 0) {
            printProgress(i, snp_count-1);
        }
        //printProgress(i, snp_count-1);
    }
    std::cout << std::endl;  // End the progress bar

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

    // Open a SAMtools process to calculate cumulative read depth and position
    // counts (non-zero depths only) for a single chromosome

    // Run the entire chromosome
        snprintf(cmd, BUFFER_SIZE,\
        "samtools depth -r %s %s | awk '{c++;s+=$3}END{print c, s}'",\
        chr.c_str(), input_filepath.c_str());  // Remove '-a' for debugging

    std::cout << "Running command: " << cmd << std::endl;

    // Parse the output
    fp = popen(cmd, "r");
    if (fp == NULL) {
        std::cerr << "ERROR: Could not open pipe for command: " << cmd << std::endl;
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
        } else {
            std::cerr << "ERROR: Could not parse output from command: " << cmd << std::endl;
            exit(EXIT_FAILURE);
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

SNPData CNVCaller::readSNPBAFs(std::string snp_filepath)
{
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

    // Extract B-allele frequency data from the VCF file
    std::cout << "Extracting allelic depth data from filtered SNPs..." << std::endl;
    cmd = "bcftools query -f '%CHROM,%POS,[%AD]\n' " + filtered_snp_vcf_filepath;
    FILE *fp = popen(cmd.c_str(), "r");
    if (fp == NULL)
    {
        std::cerr << "ERROR: Could not open pipe for command: " << cmd << std::endl;
        exit(1);
    }

    // Read the DP and AD values (AD for alternate allele depth)
    char line[BUFFER_SIZE];
    std::vector<int64_t> locations;
    std::vector<double> bafs;
    while (fgets(line, BUFFER_SIZE, fp) != NULL)
    {
        // Parse the line
        char *tok = strtok(line, ",");  // Tokenize the line
        int col = 0;  // Column index
        std::string chr = "";
        uint64_t pos = 0;
        int ref_ad = 0;
        int alt_ad = 0;
        std::string alt_allele = "";
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

            // Get the AD for the reference allele from column 3
            else if (col == 2)
            {
                ref_ad = atoi(tok);
            }

            // Get the AD for the non-reference allele from column 4
            else if (col == 3)
            {
                alt_ad = atoi(tok);
            }

            // Move to the next token
            tok = strtok(NULL, ",");
            col++;
        }

        // Calculate the BAF
        double baf = (double) alt_ad / (double) (ref_ad + alt_ad);

        // Store the position, alternate allele, and BAF
        locations.push_back(pos);
        bafs.push_back(baf);
    }

    // Close the pipe
    pclose(fp);

    // Exit if no SNPs were found
    if (locations.size() == 0)
    {
        std::cerr << "ERROR: No SNPs found in the VCF file: " << snp_filepath << std::endl;
        exit(1);
    }

    // Create the SNP data struct
    SNPData snp_data;
    snp_data.locations = locations;
    snp_data.bafs      = bafs;

    return snp_data;
}

std::vector<double> CNVCaller::getSNPPopulationFrequencies(std::vector<int64_t> snp_locations)
{
    // Get the PFB filepath
    std::string pfb_filepath = this->input_data->getPFBFilepath();

    // Open the PFB file
    FILE *fp = fopen(pfb_filepath.c_str(), "r");
    if (fp == NULL)
    {
        std::cerr << "ERROR: Could not open PFB file: " << pfb_filepath << std::endl;
        exit(1);
    }

    // Read the AF values (AF for population frequency) into a map
    std::map<std::string, std::map<int, double>> pfb_map;  // Map of population frequencies by SNP position (chr -> pos -> pfb)
    char line[BUFFER_SIZE];
    int pfb_count   = 0;
    int pfb_min_hit = 0;
    int pfb_max_hit = 0;
    while (fgets(line, BUFFER_SIZE, fp) != NULL)
    {
        // Parse the tab-delimited line
        char *tok = strtok(line, "\t");  // Tokenize the line
        int col = 0;  // Column index
        std::string chr = "";
        uint64_t pos = 0;
        double pfb = 0;
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

            // Get the PFB from column 3
            else if (col == 2)
            {
                pfb = atof(tok);
            }

            tok = strtok(NULL, "\t");
            col++;
        }

        // Store the PFB value in the map, and fix within the range
        if (pfb < MIN_PFB)
        {
            pfb_map[chr][pos] = MIN_PFB;
            pfb_min_hit++;
        } else if (pfb > MAX_PFB) {
            pfb_map[chr][pos] = MAX_PFB;
            pfb_max_hit++;
        } else {
            pfb_map[chr][pos] = pfb;

            // Print values for debugging
            //std::cout << "PFB at " << chr << ":" << pos << " -> " << std::fixed << std::setprecision(6) << pfb << std::endl;
        }
        pfb_count++;
    }

    // Close the file
    fclose(fp);

    // Log the percentage of PFB values that were fixed
    std::cout << "PFB values fixed: " << (double) (pfb_min_hit + pfb_max_hit) / (double) pfb_count * 100 << "%" << std::endl;
    std::cout << "PFB min hit: " << (double) pfb_min_hit / (double) pfb_count * 100 << "%" << std::endl;
    std::cout << "PFB max hit: " << (double) pfb_max_hit / (double) pfb_count * 100 << "%" << std::endl;

    // Determine whether the chromosome is in chr notation (e.g. chr1) or not (e.g. 1)
    bool chr_notation = false;
    std::string first_chr = pfb_map.begin()->first;
    if (first_chr.find("chr") != std::string::npos)
    {
        chr_notation = true;
    }

    // Create a vector of population frequencies for each SNP
    std::cout << "Getting SNPs population frequencies..." << std::endl;
    std::vector<double> snp_pfb;

    // Get the chromosome
    std::string chr = this->input_data->getRegionChr();

    // Make sure the chromosome follows the same format as the PFB file
    if (chr_notation)
    {
        // Add "chr" to the chromosome if it is not already there
        if (chr.find("chr") == std::string::npos)
        {
            chr = "chr" + chr;
        }
    } else {
        // Remove "chr" from the chromosome if it is there
        if (chr.find("chr") != std::string::npos)
        {
            chr = chr.substr(3);
        }
    }

    // Loop through each SNP position and get the population frequency
    int snp_count = (int) snp_locations.size();
    int found_count = 0;
    for (int i = 0; i < snp_count; i++)
    {
        // Get the SNP position
        int pos = snp_locations[i];

        // Get the population frequency from the map
        auto it = pfb_map[chr].find(pos);
        if (it != pfb_map[chr].end())
        {
            // Store the population frequency
            snp_pfb.push_back(it->second);
            found_count++;
        }
        else
        {
            // Use 0.5 as the population frequency if it is not found
            snp_pfb.push_back(0.5);
        }
    }

    // Print the percentage of SNPs with population frequencies
    std::cout << "Found population frequencies for " << found_count << " of " << snp_count << " SNPs (" << (double) found_count / (double) snp_count * 100 << "%)" << std::endl;

    return snp_pfb;
}

void CNVCaller::saveToTSV(SNPData& snp_data, std::string filepath)
{
    // Open the TSV file for writing
    std::ofstream tsv_file(filepath);

    // Write the header
    tsv_file << "chromosome\tposition\tb_allele_freq\tlog2_ratio\tcnv_state\tpopulation_freq" << std::endl;

    // Write the data
    std::string chr = this->input_data->getRegionChr();
    int snp_count = (int) snp_data.locations.size();
    for (int i = 0; i < snp_count; i++)
    {
        // Get the SNP data
        int64_t pos        = snp_data.locations[i];
        double  pfb        = snp_data.pfbs[i];
        double  baf        = snp_data.bafs[i];
        double  log2_ratio = snp_data.log2_ratios[i];
        int     cn_state   = snp_data.state_sequence[i];

        // Write the TSV line (chrom, pos, baf, lrr, state)
        tsv_file << \
            chr          << "\t" << \
            pos          << "\t" << \
            baf          << "\t" << \
            log2_ratio   << "\t" << \
            cn_state     << "\t" << \
            pfb          << \
        std::endl;
    }

    // Close the file
    tsv_file.close();
}

void CNVCaller::saveToBED(SNPData& snp_data, std::string filepath)
{
    // Each of the 6 copy number states corresponds to the following:
    // 1: 0/0 (Two copy loss)
    // 2: 1/0 (One copy loss)
    // 3: 1/1 (Normal)
    // 4: 1/1 (Copy neutral LOH)
    // 5: 2/1 (One copy gain)
    // 6: 2/2 (Two copy gain)

    // Create a map of descriptions and RGB colors for each CNV state
    // Descriptions:
    std::map<int, std::string> description_map = {
        {1, "DEL"},  // DEL
        {2, "DEL"},  // DEL
        {3, "NEUT"},  // NEUTRAL
        {4, "NEUT"},  // NEUTRAL
        {5, "DUP"},  // DUP
        {6, "DUP"}  // DUP
    };

    // RGB colors:
    std::map<int, std::string> color_map = {
        {1, "255,0,0"},  // DEL
        {2, "255,0,0"},  // DEL
        {3, "0,0,0"},  // NEUTRAL
        {4, "0,0,0"},  // NEUTRAL
        {5, "0,0,255"},  // DUP
        {6, "0,0,255"}  // DUP
    };

    // Open the BED file for writing
    std::ofstream bed_file(filepath);

    // Write the track line (track type, name, description)
    bed_file << "track name=\"CNV States\" description=\"CNV state predictions from ContextSV\" itemRgb=\"On\" visibility=\"2\" useScore=\"0\"" << std::endl;

    // Write the data
    std::string chr = this->input_data->getRegionChr();
    int snp_count = (int) snp_data.locations.size();
    for (int i = 0; i < snp_count; i++)
    {
        // Get the SNP position and state
        int64_t pos        = snp_data.locations[i];
        int     state   = snp_data.state_sequence[i];

        // Get the state description
        std::string state_str = description_map[state];

        // Get the RGB color
        std::string rgb_color = color_map[state];

        // Write the BED line (chrom, start, end, name (state), score, strand, thickStart, thickEnd, itemRgb)
        bed_file         << \
            chr          << "\t" << \
            pos          << "\t" << \
            pos          << "\t" << \
            state_str    << "\t" << \
            "0"          << "\t" << \
            "."          << "\t" << \
            pos          << "\t" << \
            pos          << "\t" << \
            rgb_color    << \
        std::endl;
    }

    // Close the file
    bed_file.close();
}
