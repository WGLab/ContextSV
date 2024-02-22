
#include "cnv_caller.h"

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
#include <numeric>  // std::iota
#include <thread>
#include <future>
#include <string>

#include "utils.h"
#include "sv_data.h"
#include "sv_types.h"

#define MIN_PFB 0.01
#define MAX_PFB 0.99
/// @endcond

using namespace sv_types;

void CNVCaller::runCopyNumberPrediction(std::string chr, SVData& sv_calls, SNPInfo& snp_info, SNPData& snp_data, CHMM hmm, int window_size)
{
    // Get the chromosome SV candidates
    std::map<SVCandidate, SVInfo>& sv_candidates = sv_calls.getChromosomeSVs(chr);
    int sv_count = (int) sv_candidates.size();
    printMessage("Total SV count: " + std::to_string(sv_count));

    // If there are no SV candidates, then return
    if (sv_count == 0)
    {
        std::cout << "No SV candidates found for chromosome " << chr << std::endl;
        return;
    }

    // Get the mean chromosome coverage
    double mean_chr_cov = 0;
    try
    {
        mean_chr_cov = this->input_data->getMeanChromosomeCoverage(chr);
        std::cout << "User-provided mean chromosome coverage for " << chr << ": " << mean_chr_cov << std::endl;
    }
    catch(const std::out_of_range& e)
    {
        // No user-provided mean chromosome coverage
        std::cout << "Calculating mean chromosome coverage for " << chr << "..." << std::endl;
        mean_chr_cov = calculateMeanChromosomeCoverage(chr);
        std::cout << "Mean chromosome coverage for " << chr << ": " << mean_chr_cov << std::endl;
    }
    
    // Run copy number prediction for the SV candidates
    // Loop through each SV candidate and predict the copy number state
    printMessage("Predicting copy number states for chromosome " + chr + "...");
    int current_sv = 0;
    for (auto const& sv_call : sv_candidates)
    {
        current_sv++;
        printMessage("Predicting SV " + std::to_string(current_sv) + " of " + std::to_string(sv_count) + "...");
        
        // Get the SV candidate
        const SVCandidate& candidate = sv_call.first;

        // Get the start and end positions of the SV call
        int64_t start_pos = std::get<0>(candidate);
        int64_t end_pos = std::get<1>(candidate);

        // Get the read depths for the entire region
        std::unordered_map<uint64_t, int> pos_depth_map;
        calculateDepthsForSNPRegion(chr, start_pos, end_pos, pos_depth_map);

        // Continue if there are no positions in the region
        if (pos_depth_map.size() == 0)
        {
            std::cout << "No reads found in SV region " << chr << ":" << start_pos << "-" << end_pos << std::endl;
            continue;
        }

        // Loop through the SV region, calculate the log2 ratios, and run the
        // Viterbi algorithm
        //std::cout << "Calculating log2 ratios for SV region " << chr << ":" << start_pos << "-" << end_pos << "..." << std::endl;
        std::vector<double> log2_cov;
        std::vector<int> state_sequence;
        std::vector<int64_t> pos;
        std::vector<double> baf;
        std::vector<double> pfb;
        std::vector<double> pfbs;
        std::vector<bool> is_snp;
        bool snps_found = false;
        int sv_length = end_pos - start_pos;

        // Estimate the number of positions in the SV region
        int interval_count = std::ceil((double) sv_length / (double) window_size);
        int current_interval = 0;

        // Loop as a sliding non-overlapping window across the SV region
        // Note: SV coordinates are 1-based
        for (int64_t i = start_pos; i <= end_pos; i += window_size)
        {
            // Run a sliding non-overlapping window of size window_size across
            // the SV region and calculate the log2 ratio for each window
            int64_t window_start = i;
            int64_t window_end = std::min(i + window_size - 1, end_pos);
            //std::cout << "Calculating log2 ratio for window " << chr << ":" << window_start << "-" << window_end << "..." << std::endl;

            // Use the position and depth map to calculate the log2 ratio
            double cum_depth = 0;
            int pos_count = 0;
            for (int64_t j = window_start; j <= window_end; j++)
            {
                // Check if the position is in the map
                auto it = pos_depth_map.find(j);
                if (it == pos_depth_map.end())
                {
                    std::cout << "Position " << j << " not found in depth map" << std::endl;
                    continue;
                }

                // Get the depth for the position
                int depth = pos_depth_map[j];

                // Update the position count and cumulative depth
                pos_count++;
                cum_depth += depth;
            }

            // Continue if there are no positions in the window
            if (pos_count == 0)
            {
                std::cout << "No positions found in window " << chr << ":" << window_start << "-" << window_end << std::endl;
                continue;
            }

            // Calculate the window coverage log2 ratio
            double window_mean_cov = (double) cum_depth / (double) pos_count;

            // Calculate the log2 ratio for the window
            double window_log2_ratio = log2(window_mean_cov / mean_chr_cov);

            // If result is infinite, print calculation details
            if (std::isinf(window_log2_ratio))
            {
                std::cout << "[DEBUG] Window mean coverage: " << window_mean_cov << std::endl;
                std::cout << "[DEBUG] Mean chromosome coverage: " << mean_chr_cov << std::endl;
                std::cout << "[DEBUG] Log2 ratio: " << window_log2_ratio << std::endl;
            }

            // Get the SNP info for the window
            //printMessage("Querying SNPs for window " + chr + ":" + std::to_string(window_start) + "-" + std::to_string(window_end) + "...");
            std::tuple<std::vector<int64_t>, std::vector<double>, std::vector<double>> window_snps = snp_info.querySNPs(chr, window_start, window_end);
            //printMessage("Finished querying SNPs for window " + chr + ":" + std::to_string(window_start) + "-" + std::to_string(window_end) + "...");
            std::vector<int64_t>& window_pos = std::get<0>(window_snps);  // SNP positions
            std::vector<double>& window_bafs = std::get<1>(window_snps);  // B-allele frequencies
            std::vector<double>& window_pfbs = std::get<2>(window_snps);  // Population frequencies of the B allele

            // Create the boolean vector for SNPs
            std::vector<bool> window_is_snp(window_pos.size(), true);

            // If there are no SNPs in the window, then use the default BAF and
            // PFB values, and the coverage log2 ratio
            if (window_pos.size() == 0)
            {
                // Use the window center position
                pos.push_back((window_start + window_end) / 2);
                baf.push_back(0.5);
                pfb.push_back(MIN_PFB);
                log2_cov.push_back(window_log2_ratio);
                is_snp.push_back(false);

            // If there are SNPs in the window, then use the BAF and PFB values
            } else {
                snps_found = true;

                // Add the BAFs and PFBs to the vectors
                pos.insert(pos.end(), window_pos.begin(), window_pos.end());
                baf.insert(baf.end(), window_bafs.begin(), window_bafs.end());
                pfb.insert(pfb.end(), window_pfbs.begin(), window_pfbs.end());
                log2_cov.insert(log2_cov.end(), window_pos.size(), window_log2_ratio);
                is_snp.insert(is_snp.end(), window_is_snp.begin(), window_is_snp.end());
            }

            // Update the position counter
            pos_count;

            // Print progress
            current_interval++;
            //std::cout << "Progress: " << std::to_string(current_interval) << " of " << std::to_string(interval_count) << " intervals" << std::endl;
        }

        // Run the Viterbi algorithm
        //int interval_count = (int) log2_cov.size();
        int data_count = (int) log2_cov.size();
        double *lrr_ptr = log2_cov.data();
        double *baf_ptr = baf.data();
        double *pfb_ptr = pfb.data();
        int *snpdist = NULL;
        double *logprob = NULL;
        state_sequence = testVit_CHMM(hmm, data_count, lrr_ptr, baf_ptr, pfb_ptr, snpdist, logprob);
        
        // Find the most common CNV state (1-6, 0-based index to 1-based index)
        std::vector<int> state_counts(6, 0);
        for (int i = 0; i < data_count; i++)
        {
            int state = state_sequence[i];
            state_counts[state - 1]++;
        }

        // Find the most common CNV state (1-6, 0-based index to 1-based index)
        int max_state = std::distance(state_counts.begin(), std::max_element(state_counts.begin(), state_counts.end())) + 1;

        // Update the SV calls with the CNV type and genotype
        int cnv_type = cnv_type_map[max_state];
        std::string genotype = cnv_genotype_map[max_state];
        std::string data_type = "Log2CNV";

        // Update the CNV type
        if (snps_found)
        {
            data_type = "SNPCNV";
        }
        
        // Update the SV type if it is not unknown
        if (cnv_type != sv_types::UNKNOWN)
        {
            //sv_info.sv_type = cnv_type;
            sv_calls.updateSVType(chr, candidate, cnv_type, data_type);

            std::cout << "[TEST] Found SV type: " << sv_types::SVTypeString[cnv_type] << std::endl;
        }

        // Update the SV genotype
        sv_calls.updateGenotype(chr, candidate, genotype);

        // Update the SNP data
        snp_data.pos.insert(snp_data.pos.end(), pos.begin(), pos.end());
        snp_data.baf.insert(snp_data.baf.end(), baf.begin(), baf.end());
        snp_data.log2_cov.insert(snp_data.log2_cov.end(), log2_cov.begin(), log2_cov.end());
        snp_data.state_sequence.insert(snp_data.state_sequence.end(), state_sequence.begin(), state_sequence.end());
        snp_data.pfb.insert(snp_data.pfb.end(), pfb.begin(), pfb.end());
        snp_data.is_snp.insert(snp_data.is_snp.end(), is_snp.begin(), is_snp.end());
    }

    std::cout << "Finished predicting copy number states for chromosome " << chr << std::endl;
    // std::cout << "Saved copy number data to " << output_tsv << std::endl;

    // Print the first three positions and BAFs
    std::cout << "First three positions and BAFs:" << std::endl;
    for (int i = 0; i < 3; i++)
    {
        std::cout << snp_data.pos[i] << "\t" << snp_data.baf[i] << std::endl;
    }
}


CNVCaller::CNVCaller(InputData &input_data)
{
    this->input_data = &input_data;
}

void CNVCaller::run(SVData& sv_calls)
{
    // Predict copy number states at SNP positions using a hidden Markov model
    // Get the region data
    bool whole_genome = this->input_data->getWholeGenome();
    std::vector<std::string> chromosomes;
    if (whole_genome)
    {
        chromosomes = this->input_data->getRefGenomeChromosomes();
    } else {
        chromosomes.push_back(this->input_data->getRegionChr());
    }

    // Read SNP positions and B-allele frequency values from the VCF file
    SNPInfo snp_info;
    std::cout << "Reading SNP allele frequencies from VCF file..." << std::endl;
    std::string snp_filepath = this->input_data->getSNPFilepath();
    readSNPAlleleFrequencies(snp_filepath, snp_info, whole_genome);

    // Get the population frequencies for each SNP
    std::cout << "Obtaining SNP population frequencies..." << std::endl;
    getSNPPopulationFrequencies(snp_info);

    // Read the HMM from file
    std::string hmm_filepath = this->input_data->getHMMFilepath();
    std::cout << "Reading HMM from file: " << hmm_filepath << std::endl;
    CHMM hmm = ReadCHMM(hmm_filepath.c_str());

    // Get the number of threads
    int num_threads = this->input_data->getThreadCount();

    // Get the window size for HMM observations
    int window_size = this->input_data->getWindowSize();

    // Loop over SV call chromosomes
    std::cout << "Predicting copy number states for SV candidates..." << std::endl;
    std::set<std::string> sv_chromosomes = sv_calls.getChromosomes();
    int sv_no_snp_count = 0;
    SNPData snp_data;
    for (auto const& chr : sv_chromosomes)
    {
        // Get the SV candidates for the chromosome if any
        std::cout << "Getting chromosome SVs for " << chr << "...";

        // Run the copy number prediction for the entire chromosome
        std::cout << "Running copy number prediction for chromosome " << chr << "..." << std::endl;
        runCopyNumberPrediction(chr, sv_calls, snp_info, snp_data, hmm, window_size);

        // Print progress
        std::cout << "Finished predicting copy number states for chromosome " << chr << std::endl;
    }

    std::cout << "Finished predicting copy number states for all chromosomes" << std::endl;
    std::cout << "Found a total of " << sv_calls.totalCalls() << " SV candidates" << std::endl;

    // Save the SNP data to a TSV file
    std::string snp_output_tsv = this->input_data->getOutputDir() + "/cnv_data.tsv";
    std::cout << "Saving SNP data to " << snp_output_tsv << "..." << std::endl;
    saveToTSV(snp_data, snp_output_tsv);
    std::cout << "Saved SNP data to " << snp_output_tsv << std::endl;
}

/// Calculate the mean chromosome coverage
double CNVCaller::calculateMeanChromosomeCoverage(std::string chr)
{
    std::string input_filepath = this->input_data->getShortReadBam();

    // Open a SAMtools process to calculate cumulative read depth and position
    // counts (non-zero depths only) for a single chromosome

    // Run the entire chromosome, printing the number of positions and the
    // cumulative read depth

    // Get the number of threads
    //int num_threads = std::thread::hardware_concurrency();
    int num_threads = this->input_data->getThreadCount();

    // Split the chromosome into equal parts for each thread
    int chr_len = this->input_data->getRefGenomeChromosomeLength(chr);
    //std::cout << "Chr " << chr << " length: " << chr_len << std::endl;

    // Split the chromosome into equal parts for each thread
    int chunk_size = chr_len / num_threads;
    std::vector<std::string> region_chunks;
    for (int i = 0; i < num_threads; i++)
    {
        int start = i * chunk_size + 1;  // 1-based index
        int end = start + chunk_size;
        if (i == num_threads - 1)
        {
            end = chr_len;
        }
        //std::cout << "Chromosome chunk: " << chr << ":" << start << "-" << end << std::endl;
        region_chunks.push_back(chr + ":" + std::to_string(start) + "-" + std::to_string(end));
    }

    // Loop through each region chunk and get the mean chromosome coverage in
    // parallel
    uint64_t pos_count = 0;
    uint64_t cum_depth = 0;
    std::vector<std::future<std::tuple<uint64_t, uint64_t>>> futures;
    for (const auto& region_chunk : region_chunks)
    {
        // Create a lambda function to get the mean chromosome coverage for the
        // region chunk
        auto get_mean_chr_cov = [region_chunk, input_filepath]() -> std::tuple<uint64_t, uint64_t>
        {
            // Run samtools depth on the entire region, and print positions and
            // depths (not chromosome)
            const int cmd_size = 256;
            char cmd[cmd_size];
            snprintf(cmd, cmd_size,\
                "samtools depth -r %s %s | awk '{c++;s+=$3}END{print c, s}'",\
                region_chunk.c_str(), input_filepath.c_str());

            // Open a pipe to read the output of the command
            FILE *fp = popen(cmd, "r");
            if (fp == NULL)
            {
                std::cerr << "ERROR: Could not open pipe for command: " << cmd << std::endl;
                exit(EXIT_FAILURE);
            }

            // Parse the outputs
            uint64_t pos_count, cum_depth;
            const int line_size = 256;
            char line[line_size];
            if (fgets(line, line_size, fp) != NULL)
            {           
                if (sscanf(line, "%ld%ld", &pos_count, &cum_depth) == 2)
                {
                    // Update the position count and cumulative depth
                    pos_count += pos_count;
                    cum_depth += cum_depth;
                } else {
                    // Do nothing. This is due to the region not having
                    // any reads.
                }
            }
            pclose(fp);  // Close the process

            return std::make_tuple(pos_count, cum_depth);
        };

        // Create a future for the thread
        std::future<std::tuple<uint64_t, uint64_t>> future = std::async(std::launch::async, get_mean_chr_cov);

        // Add the future to the vector
        futures.push_back(std::move(future));
    }

    // Loop through the futures and get the results
    for (auto& future : futures)
    {
        // Wait for the future to finish
        future.wait();

        // Get the result from the future
        std::tuple<uint64_t, uint64_t> result = future.get();

        // Update the position count and cumulative depth
        pos_count += std::get<0>(result);
        cum_depth += std::get<1>(result);
    }

    // Calculate the mean chromosome coverage
    double mean_chr_cov = (double) cum_depth / (double) pos_count;

    return mean_chr_cov;
}

void CNVCaller::calculateDepthsForSNPRegion(std::string chr, int start_pos, int end_pos, std::unordered_map<uint64_t, int>& pos_depth_map)
{
    std::string input_filepath = this->input_data->getShortReadBam();

    // Get the number of threads
    //int num_threads = std::thread::hardware_concurrency();
    int num_threads = this->input_data->getThreadCount();

    // Split the region into equal parts for each thread
    int region_size = end_pos - start_pos;
    int chunk_size = region_size / num_threads;
    std::vector<std::string> region_chunks;
    for (int i = 0; i < num_threads; i++)
    {
        int start = start_pos + i * chunk_size + 1;  // 1-based index
        int end = start + chunk_size;
        if (i == num_threads - 1)
        {
            end = end_pos;
        }
        //std::cout << "SNP region chunk: " << chr << ":" << start << "-" << end << std::endl;
        region_chunks.push_back(chr + ":" + std::to_string(start) + "-" + std::to_string(end));
    }

    // Loop through each region chunk and get the mean chromosome coverage in
    // parallel
    //std::cout << "Calculating read depths for SNP region " << chr << ":" << start_pos << "-" << end_pos << "..." << std::endl;
    std::vector<std::future<std::unordered_map<uint64_t, int>>> futures;
    for (const auto& region_chunk : region_chunks)
    {
        // Create a lambda function to get the mean chromosome coverage for the
        // region chunk
        auto get_pos_depth_map = [region_chunk, input_filepath]() -> std::unordered_map<uint64_t, int>
        {
            // Run samtools depth on the entire region, and print positions and
            // depths (not chromosome)
            const int cmd_size = 256;
            char cmd[cmd_size];
            snprintf(cmd, cmd_size,\
                "samtools depth -r %s %s | awk '{print $2, $3}'",\
                region_chunk.c_str(), input_filepath.c_str());

            // Open a pipe to read the output of the command
            FILE *fp = popen(cmd, "r");
            if (fp == NULL)
            {
                std::cerr << "ERROR: Could not open pipe for command: " << cmd << std::endl;
                exit(EXIT_FAILURE);
            }

            // Create a map of positions and depths
            std::unordered_map<uint64_t, int> pos_depth_map;
            const int line_size = 1024;
            char line[line_size];
            while (fgets(line, line_size, fp) != NULL)
            {
                // Parse the line
                uint64_t pos;
                int depth;
                if (sscanf(line, "%ld%d", &pos, &depth) == 2)
                {
                    // Add the position and depth to the map
                    pos_depth_map[pos] = depth;
                } else {
                    // Do nothing. This is due to the region not having
                    // any reads.
                }
            }

            // Close the pipe
            pclose(fp);

            return pos_depth_map;
        };

        // Create a future for the thread
        std::future<std::unordered_map<uint64_t, int>> future = std::async(std::launch::async, get_pos_depth_map);

        // Add the future to the vector
        futures.push_back(std::move(future));
    }

    // Loop through the futures and get the results
    for (auto& future : futures)
    {
        // Wait for the future to finish
        future.wait();

        // Get the result from the future
        std::unordered_map<uint64_t, int> result = future.get();

        // Merge the result into the map of positions and depths
        pos_depth_map.insert(result.begin(), result.end());
    }
    //std::cout << "Finished calculating read depths for SNP region " << chr << ":" << start_pos << "-" << end_pos << std::endl;

    // return pos_depth_map;
}

void CNVCaller::readSNPAlleleFrequencies(std::string snp_filepath, SNPInfo& snp_info, bool whole_genome)
{
    // Create a VCF filepath of filtered SNPs
    std::string filtered_snp_vcf_filepath = this->input_data->getOutputDir() + "/filtered_snps.vcf";
    std::cout << "Parsing SNPs from " << snp_filepath << std::endl;

    // Check that the SNP file is sorted by running bcftools index and reading
    // the error output
    std::string index_cmd = "bcftools index " + snp_filepath + " 2>&1 | grep -i error";
    if (this->input_data->getVerbose()) {
        std::cout << "Command: " << index_cmd << std::endl;
    }

    // Open a pipe to read the output of the command
    FILE *index_fp = popen(index_cmd.c_str(), "r");
    if (index_fp == NULL)
    {
        std::cerr << "ERROR: Could not open pipe for command: " << index_cmd << std::endl;
        exit(1);
    }

    // Read the output of the command
    const int error_size = 256;
    char index_error[error_size];
    while (fgets(index_error, error_size, index_fp) != NULL)
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
        cmd = "bcftools view -v snps -i 'QUAL > 30 && DP > 10 && FILTER = \"PASS\"' " + snp_filepath + " > " + filtered_snp_vcf_filepath;
    } else {
        std::cout << "Filtering SNPs by depth, quality, and region..." << std::endl;
        std::string region = this->input_data->getRegion();

        // Update the command to sort by chromosome and position
        cmd = "bcftools view -r " + region + " -v snps -i 'QUAL > 30 && DP > 10 && FILTER = \"PASS\"' " + snp_filepath + " > " + filtered_snp_vcf_filepath;
    }

    if (this->input_data->getVerbose()) {
        std::cout << "Filtering SNPs by depth and quality..." << std::endl;
        std::cout << "Command: " << cmd << std::endl;
    }
    system(cmd.c_str());
    
    if (this->input_data->getVerbose()) {
        std::cout << "Filtered SNPs written to " << filtered_snp_vcf_filepath << std::endl;
    }

    // Extract B-allele frequency data from the VCF file and sort by chromosome
    // and position
    std::cout << "Extracting allelic depth data from filtered SNPs..." << std::endl;
    cmd = "bcftools query -f '%CHROM,%POS,[%AD]\n' " + filtered_snp_vcf_filepath + " 2>/dev/null";
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
    const int line_size = 256;
    char line[line_size];  // Line buffer
    // SNPData *snp_data;   // Pointer to the SNP data for the current chromosome
    std::vector<int64_t> locations;
    std::vector<double> bafs;
    while (fgets(line, line_size, fp) != NULL)
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

        // Remove the 'chr' prefix from the chromosome name for SNP data. All
        // SNP data in this program does not use the 'chr' prefix
        std::string chr_snp = chr;
        if (chr_snp.find("chr") != std::string::npos)
        {
            chr_snp = chr_snp.substr(3);
        }

        // Add a new location and BAF value to the chromosome's SNP data
        // (population frequency and log2 ratio will be added later)
        snp_info.insertSNPAlleleFrequency(chr_snp, pos, baf);
    }

    // Close the pipe
    pclose(fp);
}

void CNVCaller::getSNPPopulationFrequencies(SNPInfo& snp_info)
{
    // Get the number of available threads
    int num_threads = this->input_data->getThreadCount();

    // Get the chromosome names
    std::vector<std::string> chr_names = snp_info.getChromosomes();

    // // Get the number of chromosomes
    int num_chromosomes = (int) chr_names.size();
    int chr_count = 0;

    // Loop through each chromosome in list
    for (const auto& chr : chr_names)
    {
        // std::string chr = pair.first;
        // SNPData& snp_data = pair.second;

        // Get the population frequency file for the chromosome
        std::string pfb_filepath = this->input_data->getAlleleFreqFilepath(chr);

        // Check if the filepath uses the 'chr' prefix notations based on the
        // chromosome name (e.g., *.chr1.vcf.gz vs *.1.vcf.gz)
        std::string chr_gnomad = chr;  // gnomAD data may or may not have the 'chr' prefix
        std::string chr_prefix = "chr";
        bool chr_prefix_flag = false;  // Flag to indicate if the 'chr' prefix is used
        if (pfb_filepath.find(chr_prefix) == std::string::npos)
        {
            // gnomaAD does not use the 'chr' prefix
            // Remove the 'chr' prefix from the chromosome name
            if (chr_gnomad.find(chr_prefix) != std::string::npos)
            {
                chr_gnomad = chr_gnomad.substr(chr_prefix.length());
            }
        } else {
            // Add the 'chr' prefix to the chromosome name
            if (chr_gnomad.find(chr_prefix) == std::string::npos)
            {
                chr_gnomad = chr_prefix + chr;
            }
        }

        // Remove the 'chr' prefix from the chromosome name for SNP data. All
        // SNP data in this program does not use the 'chr' prefix
        std::string chr_snp = chr;
        if (chr_snp.find(chr_prefix) != std::string::npos)
        {
            chr_snp = chr_snp.substr(chr_prefix.length());
        }

        // If no population frequency file is provided, use 0.5 as the
        // population frequency for all SNPs
        if (pfb_filepath == "")
        {
            std::cout << "No population frequency file provided for chromosome " << chr << std::endl;
            continue;
        }

        std::cout << "Reading population frequencies for chromosome " << chr << " from " << pfb_filepath << std::endl;

        // Get the start and end SNP positions for the chromosome (1-based
        // index)
        std::pair<int64_t, int64_t> snp_range = snp_info.getSNPRange(chr);
        int64_t snp_start = snp_range.first;
        int64_t snp_end = snp_range.second;
        // int snp_count = (int) snp_data.locations.size();
        // int snp_start = snp_data.locations[0];
        // int snp_end = snp_data.locations[snp_count - 1];

        // Get the number of avaiable threads
        //int num_threads = std::thread::hardware_concurrency();
        int num_threads = this->input_data->getThreadCount();

        // Get the region size
        int64_t region_size = snp_end - snp_start;

        // Split the region into equal parts for each thread
        int region_size_per_thread = region_size / num_threads;
        std::vector<std::string> region_chunks;
        for (int i = 0; i < num_threads; i++)
        {
            int64_t start = snp_start + i * region_size_per_thread;
            int64_t end = start + region_size_per_thread;
            region_chunks.push_back(chr_gnomad + ":" + std::to_string(start) + "-" + std::to_string(end));
        }

        // Loop through each region chunk and get the population frequencies in
        // parallel
        std::unordered_map<int, double> pos_pfb_map;
        std::vector<std::thread> threads;
        
        // Vector of futures
        std::vector<std::future<std::unordered_map<int, double>>> futures;
        for (const auto& region_chunk : region_chunks)
        {
            // Create a lambda function to get the population frequencies for the
            // region chunk
            auto get_pfb = [region_chunk, pfb_filepath]() -> std::unordered_map<int, double>
            {
                // Run bcftools query to get the population frequencies for the
                // chromosome within the SNP region, assumed to be sorted by
                // position
                std::string cmd = \
                    "bcftools query -r " + region_chunk + " -f '%POS\t%AF\n' -i 'INFO/variant_type=\"snv\"' " + pfb_filepath + " 2>/dev/null";

                //std::cout << "Command: " << cmd << std::endl;

                // Open a pipe to read the output of the command
                FILE *fp = popen(cmd.c_str(), "r");
                if (fp == NULL)
                {
                    std::cerr << "ERROR: Could not open pipe for command: " << cmd << std::endl;
                    exit(1);
                }

                // Loop through the BCFTOOLS output and populate the map of population
                // frequencies
                std::unordered_map<int, double> pos_pfb_map;
                const int line_size = 256;
                char line[line_size];
                while (fgets(line, line_size, fp) != NULL)
                {
                    // Parse the line
                    int pos;
                    double pfb;
                    if (sscanf(line, "%d%lf", &pos, &pfb) == 2)
                    {
                        // Add the position and population frequency to the map
                        pos_pfb_map[pos] = pfb;
                    }
                }

                // Close the pipe
                pclose(fp);

                return pos_pfb_map;
            };

            // Create a future for the thread
            std::future<std::unordered_map<int, double>> future = std::async(std::launch::async, get_pfb);
            futures.push_back(std::move(future));
        }

        // Loop through the futures and get the results
        //std::cout << "Merging population frequencies for chromosome " << chr << "..." << std::endl;
        for (auto& future : futures)
        {
            // Wait for the future to finish
            future.wait();

            // Get the result from the future
            std::unordered_map<int, double> result = future.get();

            // Loop through the result and add to SNPInfo
            for (auto& pair : result)
            {
                int pos = pair.first;
                double pfb = pair.second;
                snp_info.insertSNPPopulationFrequency(chr_snp, pos, pfb);
            }
        }

        // Update the progress (number of chromosomes processed)
        chr_count++;
        std::cout << std::to_string(chr_count) + "/" + std::to_string(num_chromosomes) + " chromosomes processed" << std::endl;
    }
}

void CNVCaller::saveToTSV(SNPData& snp_data, std::string filepath)
{
    // Open the TSV file for writing
    std::ofstream tsv_file(filepath);

    // Write the header
    tsv_file << "chromosome\tposition\tsnp\tb_allele_freq\tlog2_ratio\tcnv_state\tpopulation_freq" << std::endl;

    // Write the data
    std::string chr = this->input_data->getRegionChr();
    int snp_count = (int) snp_data.pos.size();
    for (int i = 0; i < snp_count; i++)
    {
        // Get the SNP data
        int64_t pos        = snp_data.pos[i];
        bool    is_snp     = snp_data.is_snp[i];
        double  pfb        = snp_data.pfb[i];
        double  baf        = snp_data.baf[i];
        double  log2_ratio = snp_data.log2_cov[i];
        int     cn_state   = snp_data.state_sequence[i];

        // Write the TSV line (chrom, pos, baf, lrr, state)
        tsv_file << \
            chr          << "\t" << \
            pos          << "\t" << \
            is_snp       << "\t" << \
            baf          << "\t" << \
            log2_ratio   << "\t" << \
            cn_state     << "\t" << \
            pfb          << \
        std::endl;
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
        int snp_count = (int) snp_data.pos.size();
        for (int i = 0; i < snp_count; i++)
        {
            // Get the SNP data
            int64_t pos        = snp_data.pos[i];
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
