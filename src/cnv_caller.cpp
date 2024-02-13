
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

#include "utils.h"
#include "sv_data.h"
#include "sv_types.h"

#define MIN_PFB 0.01
#define MAX_PFB 0.99
/// @endcond

using namespace sv_types;


CNVCaller::CNVCaller(InputData& input_data)
{
    this->input_data = &input_data;
}

void CNVCaller::run(SVData& sv_calls)
{
    // Define a map of CNV genotypes by HMM predicted state.
    // Each of the 6 state predictions corresponds to a copy number state:
    // 1: 0/0 (Two copy loss: homozygous deletion, GT: 0/0)
    // 2: 1/0 (One copy loss: heterozygous deletion, GT: 0/1)
    // 3: 1/1 (Normal diploid: no copy number change, GT: 1/1)
    // 4: 1/1 (Copy neutral LOH: no copy number change, GT: 1/1)
    // 5: 2/1 (One copy gain: heterozygous duplication, GT: 1/2)
    // 6: 2/2 (Two copy gain: homozygous duplication, GT: 2/2)
    std ::map<int, std::string> cnv_genotype_map = {
        {1, "0/0"},
        {2, "0/1"},
        {3, "1/1"},
        {4, "1/1"},
        {5, "1/2"},
        {6, "2/2"}
    };

    // Define a map of CNV types by HMM predicted state.
    std ::map<int, int> cnv_type_map = {
        {1, sv_types::DEL},
        {2, sv_types::DEL},
        {3, sv_types::UNKNOWN},
        {4, sv_types::UNKNOWN},
        {5, sv_types::DUP},
        {6, sv_types::DUP}
    };

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

    // Initialize the SNP data map (chr -> SNPData)
    SNPDataMap snp_data_map;
    for (auto const& chr : chromosomes)
    {
        snp_data_map[chr] = SNPData();
    }

    // Read SNP positions and B-allele frequency values from the VCF file
    std::cout << "Reading SNP allele frequencies from VCF file..." << std::endl;
    std::string snp_filepath = this->input_data->getSNPFilepath();
    readSNPAlleleFrequencies(snp_filepath, snp_data_map, whole_genome);

    // Get the population frequencies for each SNP
    std::cout << "Obtaining SNP population frequencies..." << std::endl;
    //PFBMap pfb_map = this->input_data->getPFBMap();
    //getSNPPopulationFrequencies(pfb_map, snp_data_map);
    getSNPPopulationFrequencies(snp_data_map);

    // Calculate LRRs
    calculateLog2RatioAtSNPS(snp_data_map);

    // Read the HMM from file
    //std::string hmm_filepath = "data/wgs.hmm";
    std::string hmm_filepath = this->input_data->getHMMFilepath();
    std::cout << "Reading HMM from file: " << hmm_filepath << std::endl;
    CHMM hmm = ReadCHMM(hmm_filepath.c_str());

    // Loop through each chromosome and run the Viterbi algorithm
    for (auto const& chr : chromosomes)
    {
        // Get the data vectors for the chromosome
        SNPData& snp_data = snp_data_map[chr];
        int snp_count = (int) snp_data.locations.size();

        // Run the Viterbi algorithm if there are SNPs
        if (snp_count > 0)
        {
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
        }
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

    // Loop through each SV call and check for overlap with CNV calls from the
    // SNP data
    std::cout << "Checking for overlap with CNV calls..." << std::endl;
    std::unordered_map<std::string, std::vector<SVCandidate>> svs_no_snps;

    // Loop over SV call chromosomes
    std::set<std::string> sv_chromosomes = sv_calls.getChromosomes();
    int sv_no_snp_count = 0;
    for (auto const& chr : sv_chromosomes)
    {
        // Get the SV candidates for the chromosome if any
        std::map<SVCandidate, SVInfo> sv_candidates = sv_calls.getChromosomeSVs(chr);

        // Iterate over the SV candidates for the chromosome
        for (auto it = sv_candidates.begin(); it != sv_candidates.end(); it++)
        {
            SVCandidate candidate = it->first;

            // Get the SV coordinates for the candidate
            int start_pos = std::get<0>(candidate);
            int end_pos = std::get<1>(candidate);

            // Get the SNP data for the chromosome
            SNPData& snp_data = snp_data_map[chr];

            // Check for overlap with CNV calls. Since the SNP data is sorted by
            // position, we can use a binary search to find the SNP positions that
            // overlap with the SV call
            int snp_count = (int) snp_data.locations.size();
            std::vector<int64_t>::iterator start_it = std::lower_bound(snp_data.locations.begin(), snp_data.locations.end(), start_pos);
            std::vector<int64_t>::iterator end_it = std::lower_bound(snp_data.locations.begin(), snp_data.locations.end(), end_pos);

            // Get the indices of the SNP positions that overlap with the SV call
            int start_idx = start_it - snp_data.locations.begin();
            int end_idx = end_it - snp_data.locations.begin();

            // If no SNPs overlap with the SV call, then add the SV call to the list
            if (start_idx == snp_count || end_idx == 0)
            {
                // Add the SV call to the list of SVs with no overlap
                svs_no_snps[chr].push_back(candidate);
                sv_no_snp_count++;
                continue;
            }

            // Loop through the SNP positions that overlap with the SV call and get
            // the CNV state
            std::vector<int> state_counts(6, 0);
            int dup_count = 0;
            int del_count = 0;
            int no_call_count = 0;
            int total_count = 0;
            for (int i = start_idx; i < end_idx; i++)
            {
                // Get the SNP position and CNV state
                int64_t pos = snp_data.locations[i];
                int state = snp_data.state_sequence[i];

                // Update the state counts (Note: State index is 0-5 instead of 1-6)
                state_counts[state - 1]++;

                // Update the SV type counts
                if (state == 5 || state == 6)
                {
                    dup_count++;
                } else if (state == 1 || state == 2)
                {
                    del_count++;
                } else {
                    no_call_count++;
                }
                total_count++;
            }

            // Find the most common CNV state (1-6, 0-based index to 1-based index)
            int max_state = std::distance(state_counts.begin(), std::max_element(state_counts.begin(), state_counts.end())) + 1;

            // Check if the SV region has duplication or deletion calls for at least
            // 50% of predictions
            int sv_type = UNKNOWN;
            if (total_count > 0)
            {
                if (dup_count >= total_count / 2 && dup_count > del_count)
                {
                    // Update the CNV type
                    //std::cout << "[SNP] Updating SV type to " << sv_types::SVTypeString[DUP] << " for chromosome " << chr << " SV " << start_pos << "-" << end_pos << "..." << std::endl;
                    sv_calls.updateSVType(chr, candidate, DUP, "SNPCNV");

                    // Update the genotype based on the most common CNV state
                    sv_calls.updateGenotype(chr, candidate, cnv_genotype_map[max_state]);

                } else if (del_count >= total_count / 2 && del_count > dup_count)
                {
                    // Update the CNV type
                    sv_calls.updateSVType(chr, candidate, DEL, "SNPCNV");

                    // Update the genotype based on the most common CNV state
                    //std::cout << "[SNP] Updating SV type to " << sv_types::SVTypeString[DEL] << " for chromosome " << chr << " SV " << start_pos << "-" << end_pos << "..." << std::endl;
                    sv_calls.updateGenotype(chr, candidate, cnv_genotype_map[max_state]);
                } else {
                    // No CNV call from SNP data, thus add the SV call to the list
                    // of SVs with no overlap
                    svs_no_snps[chr].push_back(candidate);
                }
            }
        }
    }

    std::cout << "Number of SVs with no overlap: " << sv_no_snp_count << std::endl;
    std::cout << "Number of SVs with overlap: " << sv_calls.totalCalls() - svs_no_snps.size() << std::endl;

    // Loop through each chromosome and calculate the log2 ratios in parallel
    for (auto const& chr : sv_chromosomes)
    {
        // Get the chromosome SVs, if any

        // Check if the chromosome has SVs with no overlap
        std::vector<SVCandidate> chr_svs_no_snps;
        if (svs_no_snps.find(chr) != svs_no_snps.end())
        {
            chr_svs_no_snps = svs_no_snps[chr];
        } else {
            continue;
        }

        // Get the chromosome coverage
        double mean_chr_cov = this->chr_mean_cov[chr];

        // Create a vector to store the log2 ratios
        std::vector<double> sv_log2_ratios(chr_svs_no_snps.size(), 0);

        // SV candidates are already sorted by position, so get the min and max
        // positions for the chromosome to obtain the read depth for the entire
        // region
        int min_pos = std::get<0>(chr_svs_no_snps[0]);
        int max_pos = std::get<1>(chr_svs_no_snps[chr_svs_no_snps.size() - 1]);

        // Get the read depth for the entire region
        std::cout << "Calculating read depth for chromosome " << chr << " SV region " << min_pos << "-" << max_pos << "..." << std::endl;
        std::unordered_map<uint64_t, int> pos_depth_map = calculateDepthsForSNPRegion(chr, min_pos, max_pos);
        std::cout << "Finished calculating read depth for chromosome " << chr << " SV region " << min_pos << "-" << max_pos << "..." << std::endl;

        // Loop through each SV and calculate the log2 ratio
        std::cout << "Calculating log2 ratios for chromosome " << chr << " SVs..." << std::endl;
        int sv_index = 0;
        for (auto const& candidate : chr_svs_no_snps)
        {
            // Get the SV coordinates
            int window_start = std::get<0>(candidate);
            int window_end = std::get<1>(candidate);

            // If the region is larger than the window size, then use the window
            // size centered at the SV position. Otherwise, use the entire region
            if (window_end - window_start > this->input_data->getWindowSize())
            {
                // Get the window start and end positions
                int center = (window_start + window_end) / 2;
                window_start = std::max(min_pos, center - this->input_data->getWindowSize() / 2);
                window_end = std::min(max_pos, center + this->input_data->getWindowSize() / 2);
                //printMessage("Large SV region of length " + std::to_string(window_end - window_start) + "...");
                //std::cout << "Large SV region of length " << end_pos - start_pos << "..." << std::endl;
            }

            // Loop through each position in the window and calculate the mean
            // read depth
            int pos_count = 0;
            int cum_depth = 0;
            for (int j = window_start; j <= window_end; j++)
            {
                // Get the depth for the position
                int depth = pos_depth_map[j];

                // Skip positions with no depth information (depth = 0)
                if (depth == 0)
                {
                    continue;
                }

                // Update the position count and cumulative depth
                pos_count++;
                cum_depth += depth;
            }

            // Calculate the mean window coverage
            double mean_window_cov = (double) cum_depth / (double) pos_count;

            // Calculate the LRR
            double l2r = log2(mean_window_cov / mean_chr_cov);

            // Store the LRR
            sv_log2_ratios[sv_index] = l2r;
            sv_index++;
        }

        // Run the Viterbi algorithm using the log2 ratios
        std::cout << "Running the Viterbi algorithm for chromosome " << chr << " SVs..." << std::endl;
        double *lrr_ptr = sv_log2_ratios.data();
        double *baf_ptr = std::vector<double>(chr_svs_no_snps.size(), 0.5).data();
        double *pfb_ptr = std::vector<double>(chr_svs_no_snps.size(), MIN_PFB).data();
        int *snpdist = NULL;
        double *logprob = NULL;

        std::cout << "SV Count: " << chr_svs_no_snps.size() << std::endl;

        // Run the Viterbi algorithm
        std::vector<int> sv_state_sequence = testVit_CHMM(hmm, chr_svs_no_snps.size(), lrr_ptr, baf_ptr, pfb_ptr, snpdist, logprob);
        std::cout << "Finished running the Viterbi algorithm for chromosome " << chr << " SVs..." << std::endl;

        // Update the SV calls with the CNV type and genotype
        std::cout << "Updating SV calls for chromosome " << chr << "..." << std::endl;
        int chr_sv_count = (int) chr_svs_no_snps.size();
        for (int i = 0; i < chr_sv_count; i++)
        {
            SVCandidate candidate = chr_svs_no_snps[i];
            int state = sv_state_sequence[i];
            int cnv_type = cnv_type_map[state];
            std::string genotype = cnv_genotype_map[state];

            // Update the SV type if not unknown
            if (cnv_type != UNKNOWN)
            {
                //std::cout << "Updating SV type to " << sv_types::SVTypeString[cnv_type] << " for chromosome " << chr << " SV " << i << "..." << std::endl;
                sv_calls.updateSVType(chr, candidate, cnv_type, "Log2CNV");
            }

            // Update the genotype
            //std::cout << "Updating genotype for chromosome " << chr << " SV " << i << "..." << std::endl;
            sv_calls.updateGenotype(chr, candidate, genotype);
        }
    }
}

void CNVCaller::calculateLog2RatioAtSNPS(SNPDataMap& snp_data_map)
{
    std::string input_filepath = this->input_data->getShortReadBam();

    // Loop through each chromosome in the map and calculate log2 ratios
    int window_size = this->input_data->getWindowSize();
    double mean_chr_cov = -1;
    std::string chr = "";
    for (auto& pair : snp_data_map)
    {
        chr = pair.first;
        SNPData& snp_data = pair.second;

        // Skip the chromosome if there are no SNPs
        if (snp_data.locations.size() == 0)
        {
            std::cerr << "WARNING: No SNPs found for chromosome " << chr << std::endl;
            continue;
        }

        // Check if there is a user-provided mean chromosome coverage
        try {
            mean_chr_cov = this->input_data->getMeanChromosomeCoverage(chr);
        } catch (const std::out_of_range& oor) {
            // No user-provided mean chromosome coverage
            mean_chr_cov = -1;
        }

        if (mean_chr_cov == -1)
        {
            // Calculate the mean chromosome coverage
            std::cout << "Calculating mean coverage for chromosome " << chr << "..." << std::endl;
            mean_chr_cov = calculateMeanChromosomeCoverage(chr);
            std::cout << "Mean coverage for chromosome " << chr << ": " << mean_chr_cov << std::endl;
        } else {
            std::cout << "Using user-provided mean coverage for chromosome " << chr << ": " << mean_chr_cov << std::endl;
        }

        // Set the mean chromosome coverage for the SNP data
        snp_data.mean_chr_cov = mean_chr_cov;

        // Set the mean chromosome coverage for the class
        this->chr_mean_cov[chr] = mean_chr_cov;

        // We will loop through each SNP and calculate the LRR for a window
        // centered at the SNP position. The window size is specified by the
        // user (default: 10 kb). To speed up the calculation, we will use the
        // SAMtools depth command to obtain read depth values for the entire
        // region and then calculate the mean read depth for each window.

        // Open a SAMtools process to get read depth values for the entire
        // region
        int snp_count = (int) snp_data.locations.size();
        int chr_len = this->input_data->getRefGenomeChromosomeLength(chr);
        uint64_t region_start = std::max(1, (int) snp_data.locations[0] - window_size / 2);

        // Have to check this since my github unit tests do not have a full ref
        // genome at this time.
        // TOOO: Add a full ref genome to the unit tests via a tarball
        uint64_t region_end;
        if (chr_len > 0)
        {
            region_end = std::min(chr_len, (int) snp_data.locations[snp_count - 1] + window_size / 2);
        } else {
            region_end = snp_data.locations[snp_count - 1] + window_size / 2;
        }
        std::unordered_map<uint64_t, int> pos_depth_map = calculateDepthsForSNPRegion(chr, region_start, region_end);

        // Loop through each SNP and calculate the LRR for a window centered at
        // the SNP position
        std::cout << "Calculating log2 ratios for chromosome " << chr << "..." << std::endl;
        std::vector<double> log2_ratios = std::vector<double>(snp_count, 0);
        for (int i=0; i < snp_count; i++)
        {
            // Get the SNP position
            // uint64_t pos = snp_data.locations[i];
            int pos = snp_data.locations[i];

            // Calculate the window start and end positions
            int window_start = pos - window_size / 2;
            int window_end = pos + window_size / 2;

            // Loop through each position in the window and calculate the mean
            // read depth
            int pos_count = 0;
            int cum_depth = 0;
            for (int j = window_start; j <= window_end; j++)
            {
                // Get the depth for the position
                int depth = pos_depth_map[j];

                // Skip positions with no depth information (depth = 0)
                if (depth == 0)
                {
                    continue;
                }

                // Update the position count and cumulative depth
                pos_count++;
                cum_depth += depth;
            }

            // Calculate the mean window coverage
            double mean_window_cov = (double) cum_depth / (double) pos_count;

            // Calculate the LRR
            double l2r = log2(mean_window_cov / mean_chr_cov);

            // Store the LRR
            log2_ratios[i] = l2r;
        }

        // Store the log2 ratios for the chromosome
        snp_data.log2_ratios = log2_ratios;

        std::cout << "Finished calculating log2 ratios for chromosome " << chr << std::endl;
    }
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
                    // Do nothing. This is likely due to the region not having
                    // any reads.
                    //std::cerr << "ERROR: Could not parse output from command: " << cmd << std::endl;
                    //std::cerr << "Line: " << line << std::endl;
                    //exit(EXIT_FAILURE);
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

std::unordered_map<uint64_t, int> CNVCaller::calculateDepthsForSNPRegion(std::string chr, int start_pos, int end_pos)
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
    std::unordered_map<uint64_t, int> pos_depth_map;
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
                    // Do nothing. This is likely due to the region not having
                    // any reads.
                    //std::cerr << "ERROR: Could not parse output from command: " << cmd << std::endl;
                    //std::cerr << "Line: " << line << std::endl;
                    //exit(EXIT_FAILURE);
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

    return pos_depth_map;
}

std::vector<double> CNVCaller::runBatchLog2Ratios(std::string chr, std::vector<SVCandidate> sv_candidates)
{
    int index = 0;
    int batch_size = (int) sv_candidates.size();
    std::vector<double> log2_ratios(batch_size, 0);
    for (auto const& candidate : sv_candidates)
    {
        // Get the SV coordinates
        int start_pos = std::get<0>(candidate);
        int end_pos = std::get<1>(candidate);

        // Get the mean chromosome coverage
        double mean_chr_cov = this->chr_mean_cov[chr];

        // If the region is larger than the window size, then use the window
        // size centered at the SV position. Otherwise, use the entire region
        if (end_pos - start_pos > this->input_data->getWindowSize())
        {
            // Get the window start and end positions
            int center = (start_pos + end_pos) / 2;
            start_pos = center - this->input_data->getWindowSize() / 2;
            end_pos = center + this->input_data->getWindowSize() / 2;
            printMessage("Large SV region of length " + std::to_string(end_pos - start_pos) + "...");
            //std::cout << "Large SV region of length " << end_pos - start_pos << "..." << std::endl;
        }

        // Calculate the log2 ratio for the SV
        printMessage("Calculating log2 ratio for SV " + chr + ":" + std::to_string(start_pos) + "-" + std::to_string(end_pos) + "...");
        //std::cout << "Calculating log2 ratio for SV " << chr << ":" << start_pos << "-" << end_pos << "..." << std::endl;
        double l2r = this->calculateWindowLogRRatio(mean_chr_cov, chr, start_pos, end_pos);
        printMessage("Completed " + std::to_string(index + 1) + " of " + std::to_string(batch_size) + " log2 ratios...");

        // Store the log2 ratio
        log2_ratios[index] = l2r;
        index++;
    }
    return log2_ratios;
}

double CNVCaller::calculateWindowLogRRatio(double mean_chr_cov, std::string chr, int window_start, int window_end)
{
    std::string input_filepath = this->input_data->getShortReadBam();

    // Open a SAMtools process to calculate cumulative read depth and position
    // counts (non-zero depths only) for the region
    const int cmd_size = 256;
    char cmd[cmd_size];
    FILE *fp;
    //std::cout << "Calculating log2 ratio for window " << chr << ":" << window_start << "-" << window_end << "..." << std::endl;
    snprintf(cmd, cmd_size,\
        "samtools depth -r %s:%d-%d %s | awk '{c++;s+=$3}END{print c, s}'",\
        chr.c_str(), window_start, window_end, input_filepath.c_str());

    fp = popen(cmd, "r");
    if (fp == NULL)
    {
        fprintf(stderr, "Failed to run command\n");
        exit(EXIT_FAILURE);
    }

    // Parse the outputs
    uint64_t pos_count, cum_depth;
    double region_lrr = -1;
    const int line_size = 256;
    char line[line_size];
    if (fgets(line, line_size, fp) != NULL)
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
    SNPData *snp_data;   // Pointer to the SNP data for the current chromosome
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

    // Sort the SNP data by position and get the indices to use for sorting
    std::vector<int> indices(snp_data->locations.size());
    std::iota(indices.begin(), indices.end(), 0);

    // Sort the indices by the SNP positions
    std::sort(indices.begin(), indices.end(), [&snp_data](int i, int j) { return snp_data->locations[i] < snp_data->locations[j]; });

    // Reorder the SNP data vectors by the sorted indices
    std::vector<int64_t> sorted_locations(snp_data->locations.size());
    std::vector<double> sorted_bafs(snp_data->bafs.size());
    int snp_count = (int) snp_data->locations.size();
    for (int i = 0; i < snp_count; i++)
    {
        sorted_locations[i] = snp_data->locations[indices[i]];
        sorted_bafs[i] = snp_data->bafs[indices[i]];
    }

    // Update the SNP data vectors
    snp_data->locations = sorted_locations;
    snp_data->bafs = sorted_bafs;
}

void CNVCaller::getSNPPopulationFrequencies(SNPDataMap& snp_data_map)
{
    // Get the number of chromosomes
    int num_chromosomes = (int) snp_data_map.size();
    int chr_count = 0;

    // Loop through each chromosome in the SNP data map and access the
    // population frequency for each SNP
    for (auto& pair : snp_data_map)
    {
        std::string chr = pair.first;
        SNPData& snp_data = pair.second;

        // Get the population frequency file for the chromosome
        std::string pfb_filepath = this->input_data->getAlleleFreqFilepath(chr);

        // Check if the filepath uses the 'chr' prefix notations based on the
        // chromosome name (e.g., *.chr1.vcf.gz vs *.1.vcf.gz)
        std::string chr_gnomad = chr;
        std::string chr_prefix = "chr";
        if (pfb_filepath.find(chr_prefix) == std::string::npos)
        {
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

        // If no population frequency file is provided, use 0.5 as the
        // population frequency for all SNPs
        if (pfb_filepath == "")
        {
            std::cout << "No population frequency file provided for chromosome " << chr << ". Using the minimum value of " << MIN_PFB << " for all SNPs." << std::endl;

            // Populate the PFB vector with the minimum value
            int snp_count = (int) snp_data.locations.size();
            snp_data.pfbs = std::vector<double>(snp_count, MIN_PFB);
            continue;
        }

        std::cout << "Reading population frequencies for chromosome " << chr << " from " << pfb_filepath << std::endl;

        // Get the start and end SNP positions for the chromosome (1-based index)
        int snp_count = (int) snp_data.locations.size();
        int snp_start = snp_data.locations[0];
        int snp_end = snp_data.locations[snp_count - 1];

        // Get the number of avaiable threads
        //int num_threads = std::thread::hardware_concurrency();
        int num_threads = this->input_data->getThreadCount();

        // Get the region size
        int region_size = snp_end - snp_start;

        // Split the region into equal parts for each thread
        int region_size_per_thread = region_size / num_threads;
        std::vector<std::string> region_chunks;
        for (int i = 0; i < num_threads; i++)
        {
            int start = snp_start + i * region_size_per_thread;
            int end = start + region_size_per_thread;
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

            // Merge the result into the map of population frequencies
            pos_pfb_map.insert(result.begin(), result.end());
        }
        //std::cout << "Finished merging population frequencies for chromosome " << chr << std::endl;
        //std::cout << "Found " << pos_pfb_map.size() << " population frequencies for chromosome " << chr << std::endl;
        
        // Loop through the SNP positions and set the population frequencies
        std::cout << "Setting SNP population frequencies for chromosome " << chr << "..." << std::endl;
        int snp_pos = 0;  // SNP position
        int found_count = 0;
        int min_fixed_pfb = 0;
        int max_fixed_pfb = 0;
        double pfb;  // Population frequency value
        std::vector<double> pfbs(snp_count, 0);
        for (int i = 0; i < snp_count; i++)
        {
            // Get the SNP position
            snp_pos = snp_data.locations[i];

            // If the position is found, add the population frequency to the SNP
            // data
            auto it = pos_pfb_map.find(snp_pos);
            if (it != pos_pfb_map.end())
            {
                pfb = it->second;
                if (pfb < MIN_PFB)
                {
                    pfb = MIN_PFB;
                    min_fixed_pfb++;
                } else if (pfb > MAX_PFB)
                {
                    pfb = MAX_PFB;
                    max_fixed_pfb++;
                }
                pfbs[i] = pfb;
                found_count++;
            }
        }

        // Store the population frequencies in the SNP data
        snp_data.pfbs = pfbs;

        // Print the percentage of SNPs with population frequencies
        std::cout << "Found " << found_count << " of " << snp_count << " SNP population frequencies ";
        std::cout << "(" << (double) found_count / (double) snp_count * 100 << "%)" << std::endl;
        std::cout << "Fixed " << min_fixed_pfb << " population frequencies below " << MIN_PFB << " and ";
        std::cout << max_fixed_pfb << " population frequencies above " << MAX_PFB << std::endl;


        // Update the progress (number of chromosomes processed)
        chr_count++;
        std::cout << std::to_string(chr_count) + "/" + std::to_string(num_chromosomes) + " chromosomes processed" << std::endl;
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
