
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

void CNVCaller::run(CNVData& cnv_data)
{
    // Predict copy number states at SNP positions.
    // Each of the 6 state predictions corresponds to a copy number state:
    // 1: 0/0 (Two copy loss)
    // 2: 1/0 (One copy loss)
    // 3: 1/1 (Normal)
    // 4: 1/1 (Copy neutral LOH)
    // 5: 2/1 (One copy gain)
    // 6: 2/2 (Two copy gain)

    // Get the region data
    bool whole_genome = this->input_data->getWholeGenome();
    std::vector<std::string> chromosomes;
    if (whole_genome)
    {
        chromosomes = this->input_data->getRefGenome().getChromosomes();
    } else {
        chromosomes.push_back(this->input_data->getRegionChr());
    }

    // Initialize the SNP data map (chr -> SNPData)
    SNPDataMap snp_data_map;
    for (auto const& chr : chromosomes)
    {
        snp_data_map[chr] = SNPData();
    }

    // Read SNP positions and B-allele frequency values from the VCF file
    std::cout << "Reading SNP positions and BAF values from the VCF file..." << std::endl;
    std::string snp_filepath = this->input_data->getSNPFilepath();
    readSNPAlleleFrequencies(snp_filepath, snp_data_map, whole_genome);
    std::cout << "Complete." << std::endl;
    // int snp_count = (int) snp_data.locations.size();
    // std::cout << "Complete. Number of SNPs: " << snp_count << std::endl;

    // Get the population frequencies for each SNP
    if (this->input_data->getAlleleFreqFilepaths() == "")
    {
        std::cout << "No SNP population frequency filepaths provided. Using the minimum PFB value of " << MIN_PFB << " for all SNPs." << std::endl;

        // Populate the PFB vector with the minimum value
        for (auto const& chr : chromosomes)
        {
            SNPData& snp_data = snp_data_map[chr];
            int snp_count = (int) snp_data.locations.size();
            snp_data.pfbs = std::vector<double>(snp_count, MIN_PFB);
        }

    } else {
        std::cout << "Using provided SNP population frequency filepaths." << std::endl;
        std::cout << "Getting population frequencies for each SNP..." << std::endl;
        PFBMap pfb_map = this->input_data->getPFBMap();
        getSNPPopulationFrequencies(pfb_map, snp_data_map);
        std::cout << "Complete." << std::endl;
    }

    // Calculate LRRs
    std::cout << "Calculating log2 ratio of coverage at SNP positions..." << std::endl;
    calculateLog2RatioAtSNPS(snp_data_map);
    std::cout << "Log2 ratios calculated." << std::endl;

    // Read the HMM from file
    //std::string hmm_filepath = "data/wgs.hmm";
    std::string hmm_filepath = this->input_data->getHMMFilepath();
    std::cout << "Reading HMM from file: " << hmm_filepath << std::endl;
    CHMM hmm = ReadCHMM(hmm_filepath.c_str());
    std::cout << "HMM read from file." << std::endl;

    // Loop through each chromosome and run the Viterbi algorithm
    for (auto const& chr : chromosomes)
    {
        // Get the data vectors for the chromosome
        SNPData& snp_data = snp_data_map[chr];
        int snp_count = (int) snp_data.locations.size();
        double *lrr_ptr = snp_data.log2_ratios.data();
        double *baf_ptr = snp_data.bafs.data();
        double *pfb_ptr = snp_data.pfbs.data();
        int *snpdist = NULL;
        double *logprob = NULL;

        // Run the Viterbi algorithm
        std::cout << "Running the Viterbi algorithm for chromosome " << chr << "..." << std::endl;
        std::vector<int> state_sequence;  // Create the output state sequence
        state_sequence = testVit_CHMM(hmm, snp_count, lrr_ptr, baf_ptr, pfb_ptr, snpdist, logprob);
        snp_data.state_sequence = state_sequence;
        std::cout << "Complete." << std::endl;
    }

    // Save a TSV of the positions, LRRs, BAFs, and state sequence
    std::cout << "Saving TSV of copy number prediction data..." << std::endl;
    std::string output_tsv = this->input_data->getOutputDir() + "/cnv_data.tsv";
    saveToTSV(snp_data_map, output_tsv);
    std::cout << "Saved to: " << output_tsv << std::endl;

    // Save a BED file of the CNV state sequence
    std::cout << "Saving BED of CNV states..." << std::endl;
    std::string output_bed = this->input_data->getOutputDir() + "/cnv_states.bed";
    saveToBED(snp_data_map, output_bed);
    std::cout << "Saved to: " << output_bed << std::endl;

    // // Get the chromosome as a C string
    // char *chr_cstr = new char[chr.length() + 1];
    // strcpy(chr_cstr, chr.c_str());

    // // Create a map of the state sequence by position
    // for (int i = 0; i < snp_count; i++)
    // {
    //     cnv_data.addCNVCall(chr_cstr, snp_data.locations[i], snp_data.state_sequence[i]);
    // }

    // Loop through each chromosome and add the CNV calls to the CNVData object
    for (auto const& chr : chromosomes)
    {
        // Get the SNP data for the chromosome
        SNPData& snp_data = snp_data_map[chr];

        // Get the SNP count
        int snp_count = (int) snp_data.locations.size();
        for (int i = 0; i < snp_count; i++)
        {
            // Get the SNP data
            int64_t pos        = snp_data.locations[i];
            // double  pfb        = snp_data.pfbs[i];
            // double  baf        = snp_data.bafs[i];
            // double  log2_ratio = snp_data.log2_ratios[i];
            int     cn_state   = snp_data.state_sequence[i];

            // Add the CNV call to the CNVData object
            cnv_data.addCNVCall(chr, pos, cn_state);
        }
    }

    // // Free the C string
    // delete[] chr_cstr;
}

void CNVCaller::calculateLog2RatioAtSNPS(SNPDataMap& snp_data_map)
{
    std::string input_filepath = this->input_data->getShortReadBam();
    // std::string chr = this->input_data->getRegionChr();

    // Check if the chromosome coverage was passed in for each chromosome
    std::map<std::string, double> chr_cov_map;
    int window_size = this->input_data->getWindowSize();

    // Loop through each chromosome in the map and calculate log2 ratios
    for (auto& pair : snp_data_map)
    {
        std::string chr = pair.first;
        SNPData& snp_data = pair.second;
        double mean_chr_cov = -1;
        if (this->input_data->getChrCov(chr, mean_chr_cov) == -1)
        {
            // Calculate the mean chromosome coverage
            std::cout <<  "\nCalculating coverage for chromosome: " << chr << std::endl;
            mean_chr_cov = calculateMeanChromosomeCoverage(chr);
            std::cout << "Mean coverage for chromosome " << chr << ": " << mean_chr_cov << std::endl;
        } else {
            std::cout << "Using user-provided mean coverage for chromosome " << chr << ": " << mean_chr_cov << std::endl;
        }

        // Add the mean chromosome coverage to the map
        chr_cov_map[chr] = mean_chr_cov;
    
        // Set the region start and end from the first and last SNPs
        std::cout << "Getting region start and end..." << std::endl;
        int region_start = snp_data.locations.front();
        int region_end = snp_data.locations.back();

        // Update the region start and end if a region was specified
        if (this->input_data->getRegionSet()) {
            region_start = std::max(region_start, this->input_data->getRegionStart());
            region_end = std::min(region_end, this->input_data->getRegionEnd());
        }

        // Loop through each SNP and calculate the LRR
        std::cout << "Calculating log2 ratios for each SNP at chromosome " << chr << "..." << std::endl;
        int snp_count = (int) snp_data.locations.size();
        for (int i = 0; i < snp_count; i++) {
            int pos = snp_data.locations[i];

            // Skip SNPs outside of the region
            if (pos < region_start || pos > region_end) {
                continue;
            }

            // Calculate window mean coverage
            int window_start = pos - (window_size / 2);
            int window_end = pos + (window_size / 2);
            double lrr = calculateWindowLogRRatio(mean_chr_cov, window_start, window_end);

            // Set the LRR value
            snp_data.log2_ratios.push_back(lrr);

            // Update the progress bar every 1000 SNPs
            if (i % 1000 == 0) {
                printProgress(i, snp_count-1);
            }
        }

        // End the progress bar
        printProgress(snp_count-1, snp_count-1);
        std::cout << std::endl;
    }
}

/// Calculate the mean chromosome coverage
double CNVCaller::calculateMeanChromosomeCoverage(std::string chr)
{
    std::string input_filepath = this->input_data->getShortReadBam();

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
    std::string input_filepath = this->input_data->getShortReadBam();

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

void CNVCaller::readSNPAlleleFrequencies(std::string snp_filepath, SNPDataMap& snp_data_map, bool whole_genome)
{
    // Create a VCF filepath of filtered SNPs
    std::string filtered_snp_vcf_filepath = this->input_data->getOutputDir() + "/filtered_snps.vcf";
    std::cout << "Parsing SNPs from " << snp_filepath << std::endl;

    // Check that the SNP file is sorted by running bcftools index
    std::cout << "Checking if the SNP file is sorted..." << std::endl;
    std::string index_cmd = "bcftools index " + snp_filepath + " 2>&1 | grep -i error";
    std::cout << "> " << index_cmd << std::endl;

    // Open a pipe to read the output of the command
    FILE *index_fp = popen(index_cmd.c_str(), "r");
    if (index_fp == NULL)
    {
        std::cerr << "ERROR: Could not open pipe for command: " << index_cmd << std::endl;
        exit(1);
    }

    // Read the output of the command
    char index_error[BUFFER_SIZE];
    while (fgets(index_error, BUFFER_SIZE, index_fp) != NULL)
    {
        std::cerr << "ERROR: " << index_error << std::endl;
        exit(1);
    }

    // Close the pipe
    pclose(index_fp);

    // Filter variants by depth and quality and SNPs only
    std::string cmd;
    if (whole_genome)
    {
        std::cout << "Filtering SNPs by depth and quality..." << std::endl;
        cmd = "bcftools view -v snps -i 'QUAL > 30 && DP > 10 && FILTER = \"PASS\"' " + snp_filepath + " > " + filtered_snp_vcf_filepath;
    } else {
        std::cout << "Filtering SNPs by depth, quality, and region..." << std::endl;
        std::string region = this->input_data->getRegion();

        // Update the command to sort by chromosome and position
        cmd = "bcftools view -r " + region + " -v snps -i 'QUAL > 30 && DP > 10 && FILTER = \"PASS\"' " + snp_filepath + " > " + filtered_snp_vcf_filepath;
    }

    std::cout << "> " << cmd << std::endl;
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

    // Read the DP and AD values (AD for alternate allele depth) and store in a
    // map by chromosome and position
    std::string chr = "";  // Current chromosome
    std::string alt_allele = "";  // Alternate allele
    uint64_t pos = 0;
    int ref_ad = 0;
    int alt_ad = 0;
    char line[BUFFER_SIZE];  // Line buffer
    SNPData *snp_data;   // Pointer to the SNP data for the current chromosome
    while (fgets(line, BUFFER_SIZE, fp) != NULL)
    {
        // Parse the line
        char *tok = strtok(line, ",");  // Tokenize the line
        int col = 0;  // Column index
        while (tok != NULL)
        {
            // Get the chromosome from column 1
            if (col == 0)
            {
                chr = tok;

                // Find the SNP data for the chromosome and skip if not found
                auto it = snp_data_map.find(chr);
                if (it == snp_data_map.end())
                {
                    // Skip the line if the chromosome is not found
                    tok = strtok(NULL, ",");
                    col++;
                    continue;
                }
                snp_data = &snp_data_map[chr];
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

        // Add a new location and BAF value to the chromosome's SNP data
        snp_data->locations.push_back(pos);
        snp_data->bafs.push_back(baf);
    }

    // Close the pipe
    pclose(fp);
}

void CNVCaller::getSNPPopulationFrequencies(PFBMap& pfb_map, SNPDataMap& snp_data_map)
{
    // Loop through each chromosome in the SNP data map and access the
    // population frequency for each SNP
    std::cout << "Getting SNPs population frequencies..." << std::endl;
    for (auto& pair : snp_data_map)
    {
        std::string chr = pair.first;
        SNPData& snp_data = pair.second;

        // Check if there are any SNPs for the chromosome
        if (snp_data.locations.size() == 0)
        {
            std::cerr << "WARNING: No SNPs found for chromosome " << chr << std::endl;
            continue;
        }

        // Check if the chromosome is in the PFB map
        auto it = pfb_map.find(chr);
        if (it == pfb_map.end())
        {
            std::cerr << "WARNING: No population frequencies found for chromosome " << chr << std::endl;
            continue;
        }

        // Get the SNP count
        int snp_count = (int) snp_data.locations.size();
        int found_count = 0;
        for (int i = 0; i < snp_count; i++)
        {
            // Get the SNP position
            int pos = snp_data.locations[i];

            // Get the population frequency from the map
            auto it = pfb_map[chr].find(pos);
            if (it != pfb_map[chr].end())
            {
                // Store the population frequency
                snp_data.pfbs.push_back(it->second);
                found_count++;
            }
            else
            {
                // Use 0.5 as the population frequency if it is not found
                snp_data.pfbs.push_back(0.5);
            }
        }

        // Print the percentage of SNPs with population frequencies
        std::cout << "Chromosome " << chr << ":" << std::endl;
        std::cout << "Found population frequencies for " << found_count << " of " << snp_count << " SNPs (" << (double) found_count / (double) snp_count * 100 << "%)" << std::endl;
    }
}

void CNVCaller::saveToTSV(SNPDataMap& snp_data_map, std::string filepath)
{
    // Open the TSV file for writing
    std::ofstream tsv_file(filepath);

    // Write the header
    tsv_file << "chromosome\tposition\tb_allele_freq\tlog2_ratio\tcnv_state\tpopulation_freq" << std::endl;

    // Loop through each SNP and write the data
    for (auto& pair : snp_data_map)
    {
        // Get the SNP data for the chromosome
        const std::string chr = pair.first;
        SNPData& snp_data = pair.second;

        // Get the SNP count
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
    }

    // Close the file
    tsv_file.close();
}

void CNVCaller::saveToBED(SNPDataMap& snp_data_map, std::string filepath)
{
    // Save the CNV state sequence to a BED file.
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

    // Loop through each SNP and write the data
    for (auto& pair : snp_data_map)
    {
        // Get the SNP data for the chromosome
        const std::string chr = pair.first;
        SNPData& snp_data = pair.second;

        // Get the SNP count
        int snp_count = (int) snp_data.locations.size();
        for (int i = 0; i < snp_count; i++)
        {
            // Get the SNP data
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
    }

    // Close the file
    bed_file.close();
}
