
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

// Function to call the Viterbi algorithm for the CHMM
std::vector<int> runViterbi(CHMM hmm, SNPData& snp_data)
{
    int data_count = (int) snp_data.pos.size();
    double *lrr_ptr = snp_data.log2_cov.data();
    double *baf_ptr = snp_data.baf.data();
    double *pfb_ptr = snp_data.pfb.data();
    int *snpdist = NULL;
    double *logprob = NULL;
    std::vector<int> state_sequence = testVit_CHMM(hmm, data_count, lrr_ptr, baf_ptr, pfb_ptr, snpdist, logprob);
    return state_sequence;
}

// Function to obtain SNP information for a region
std::pair<SNPData, bool> CNVCaller::querySNPRegion(std::string chr, int64_t start_pos, int64_t end_pos, SNPInfo& snp_info, std::unordered_map<uint64_t, int>& pos_depth_map, double mean_chr_cov)
{
    // Get the SNP information for the region
    SNPData snp_data;
    bool snps_found = false;
    int window_size = this->input_data->getWindowSize();

    // // If the window size is greater than 10kb, then print a message
    // if (end_pos - start_pos > 10000)
    // {
    //     printMessage("Querying SNPs for large region of size (kb): " + std::to_string((end_pos - start_pos + 1) / 1000));
    // }

    for (int64_t i = start_pos; i <= end_pos; i += window_size)
    {
        // Run a sliding non-overlapping window of size window_size across
        // the SV region and calculate the log2 ratio for each window
        int64_t window_start = i;
        int64_t window_end = std::min(i + window_size - 1, end_pos);

        // Get the SNP info for the window
        //printMessage("Finding SNPs...");
        this->snp_data_mtx.lock();
        std::tuple<std::vector<int64_t>, std::vector<double>, std::vector<double>> window_snps = snp_info.querySNPs(chr, window_start, window_end);
        this->snp_data_mtx.unlock();
        //std::cout << "Found " << std::get<0>(window_snps).size() << " SNPs in window " << chr << ":" << window_start << "-" << window_end << std::endl;
        std::vector<int64_t>& snp_window_pos = std::get<0>(window_snps);  // SNP positions
        std::vector<double>& snp_window_bafs = std::get<1>(window_snps);  // B-allele frequencies
        std::vector<double>& snp_window_pfbs = std::get<2>(window_snps);  // Population frequencies of the B allele

        // Loop though the SNP positions and calculate the log2 ratio for
        // the window up to the SNP, then calculate the log2 ratio centered
        // at the SNP, and finally calculate the log2 ratio for the window
        // after the SNP, and continue until the end of the window
        std::vector<double> window_log2_ratios;
        int snp_count = (int) snp_window_pos.size();

        // If there are no SNPs in the window, then use the default BAF and
        // PFB values, and the coverage log2 ratio
        if (snp_count == 0)
        {
            // Calculate the log2 ratio for the window
            //printMessage("Calculating log2 ratio, length (kb): " + std::to_string((window_end - window_start + 1) / 1000));
            double window_log2_ratio = calculateLog2Ratio(window_start, window_end, pos_depth_map, mean_chr_cov);

            // Use the window center position
            //printMessage("Updating SNP data...");
            double pfb_default = 0.5;
            double baf_default = 0.5;
            this->updateSNPData(snp_data, (window_start + window_end) / 2, pfb_default, baf_default, window_log2_ratio, false);

        } else {
            snps_found = true;

            // Loop through the SNPs and calculate the log2 ratios
            int64_t bin_start = window_start;
            int64_t bin_end = 0;
            for (int j = 0; j < snp_count; j++)
            {
                // Get the SNP position
                int64_t snp_pos = snp_window_pos[j];

                // SNP bin starts at 1/2 the distance between the previous SNP
                // and the current SNP, and ends at 1/2 the distance between
                // the current SNP and the next SNP. For the first SNP, the
                // bin starts at the window start and ends at 1/2 the distance
                // between the first SNP and the next SNP, and for the last
                // SNP, the bin starts at 1/2 the distance between the previous
                // SNP and the last SNP and ends at the window end.
                bin_end = snp_pos + (j == snp_count-1 ? (window_end - snp_pos) / 2 : (snp_window_pos[j+1] - snp_pos) / 2);

                // Calculate the log2 ratio for the SNP bin
                //printMessage("Calculating log2 ratio...");
                double bin_cov = calculateLog2Ratio(bin_start, bin_end, pos_depth_map, mean_chr_cov);

                //printMessage("Updating SNP data...");
                this->updateSNPData(snp_data, (bin_start + bin_end) / 2, snp_window_pfbs[j], snp_window_bafs[j], bin_cov, true);

                // Update the previous bin start
                bin_start = bin_end + 1;
            }
        }
    }

    return std::make_pair(snp_data, snps_found);
}
SNPData CNVCaller::runCopyNumberPrediction(std::string chr, std::map<SVCandidate, SVInfo>& sv_candidates, SNPInfo& snp_info, CHMM hmm, int window_size)
{
    SNPData snp_data;

    // Get the chromosome SV candidates
    //std::map<SVCandidate, SVInfo>& sv_candidates = sv_calls.getChromosomeSVs(chr);
    int sv_count = (int) sv_candidates.size();
    printMessage("Total SV count: " + std::to_string(sv_count));

    // If there are no SV candidates, then return
    if (sv_count == 0)
    {
        printMessage("No SV candidates found for chromosome " + chr);
        return snp_data;
    }

    // Get the first and last positions of the chromosome SVs
    int64_t first_pos = std::get<0>(sv_candidates.begin()->first);  // Start position is element 0
    int64_t last_pos = std::get<1>(sv_candidates.rbegin()->first);  // End position is element 1

    // Generate the read depth map for the entire SV region
    std::unordered_map<uint64_t, int> pos_depth_map;
    printMessage("Calculating read depths for SV region " + chr + ":" + std::to_string(first_pos) + "-" + std::to_string(last_pos) + "...");
    calculateDepthsForSNPRegion(chr, first_pos, last_pos, pos_depth_map);

    // Get the mean chromosome coverage
    double mean_chr_cov = 0;
    try
    {
        mean_chr_cov = this->input_data->getMeanChromosomeCoverage(chr);
        printMessage("User-provided mean chromosome coverage for " + chr + ": " + std::to_string(mean_chr_cov));
    }
    catch(const std::out_of_range& e)
    {
        // No user-provided mean chromosome coverage
        printMessage("Calculating mean chromosome coverage for " + chr + "...");
        mean_chr_cov = calculateMeanChromosomeCoverage(chr);
        printMessage("Mean chromosome coverage for " + chr + ": " + std::to_string(mean_chr_cov));
    }

    // Get whether to extend the SNP CNV regions around the SV breakpoints
    bool extend_cnv_regions = this->input_data->getExtendCNVRegions();
    
    // Run copy number prediction for the SV candidates
    // Loop through each SV candidate and predict the copy number state
    printMessage("Predicting copy number states for chromosome " + chr + "...");
    int current_sv = 0;

    // Create a map with counts for each CNV type
    std::map<int, int> cnv_type_counts;
    for (int i = 0; i < 6; i++)
    {
        cnv_type_counts[i] = 0;
    }

    int seqprint_count = 1;
    for (auto const& sv_call : sv_candidates)
    {
        current_sv++;
        // Print every 10,000th SV
        if (current_sv % 10000 == 0)
        {
            printMessage("Predicting SV " + std::to_string(current_sv) + " of " + std::to_string(sv_count) + "...");
        }
//        printMessage("Predicting SV " + std::to_string(current_sv) + " of " + std::to_string(sv_count) + "...");
        
        // Get the SV candidate
        const SVCandidate& candidate = sv_call.first;

        // Get the start and end positions of the SV call
        int64_t start_pos = std::get<0>(candidate);
        int64_t end_pos = std::get<1>(candidate);

        // Loop through the SV region, calculate the log2 ratios, and run the
        // Viterbi algorithm to predict the copy number states
        //printMessage("Querying SNPs...");
        std::pair<SNPData, bool> snp_call = this->querySNPRegion(chr, start_pos, end_pos, snp_info, pos_depth_map, mean_chr_cov);
        SNPData& sv_snps = snp_call.first;
        bool snps_found = snp_call.second;

        // Run the Viterbi algorithm
        //printMessage("Running Viterbi algorithm...");
        std::vector<int> state_sequence = runViterbi(hmm, sv_snps);

        // // Print the first state sequence
        // if (current_sv == 1)
        // {
        //     std::string state_seq_str = "State sequence: ";
        //     for (int i = 0; i < (int) state_sequence.size(); i++)
        //     {
        //         state_seq_str += std::to_string(state_sequence[i]) + " ";
        //     }
        //     printMessage("[TEST] First state sequence: " + state_seq_str);
        // }

        // // Find the most common CNV state (1-6, 0-based index to 1-based index)
        // //printMessage("Finding the most common CNV state...");
        // std::vector<int> state_counts(6, 0);
        // int data_count = (int) state_sequence.size();
        // for (int i = 0; i < data_count; i++)
        // {
        //     int state = state_sequence[i];
        //     state_counts[state - 1]++;
        // }

        // // Find the most common CNV state (1-6, 0-based index to 1-based index)
        // // int max_state = std::distance(state_counts.begin(),
        // // std::max_element(state_counts.begin(), state_counts.end())) + 1;
        
        // Determine if there is a majority state and if it is greater than 50%
        int max_state = 0;
        int max_count = 0;
        for (int i = 0; i < 6; i++)
        {
            int state_count = std::count(state_sequence.begin(), state_sequence.end(), i+1);
            if (state_count > max_count)
            {
                max_state = i+1;
                max_count = state_count;
            }
        }

        // If there is no majority state, then set the state to unknown
        int state_count = (int) state_sequence.size();
        if (max_count < state_count / 2)
        {
            max_state = 0;
        }

        // Update the SV calls with the CNV type and genotype
        int cnv_type = cnv_type_map[max_state];
        std::string genotype = cnv_genotype_map[max_state];

        // Determine the SV calling method used to call the SV
        // (SNPCNV=SNP-based, Log2CNV=coverage-based)
        std::string data_type;
        if (snps_found)
        {
            data_type = "SNPCNV";
        } else {
            data_type = "Log2CNV";
        }

        // [TEMPORARY: Commented out for testing purposes]        
        // Update the SV type if it is not unknown
        if (cnv_type != sv_types::UNKNOWN)
        {
            this->updateSVType(sv_candidates, sv_call.first, cnv_type, data_type);
            // sv_candidates[sv_call.first].sv_type = cnv_type;

            // // [TEST] If deletion, print the state sequence
            // if (cnv_type == sv_types::DEL && seqprint_count < 10)
            // {
            //     std::string state_seq_str = "DELSEQ: ";
            //     for (int i = 0; i < (int) state_sequence.size(); i++)
            //     {
            //         state_seq_str += std::to_string(state_sequence[i]) + " ";
            //     }
            //     printMessage(state_seq_str);
            //     seqprint_count++;
            // }

            // Update the SV type counts
            cnv_type_counts[cnv_type]++;
        }

        // Update the SV genotype
        // sv_calls.updateGenotype(chr, candidate, genotype);
        this->updateSVGenotype(sv_candidates, sv_call.first, genotype);
        // sv_candidates[sv_call.first].genotype = genotype;

        // Preceding SNPs:
        int64_t sv_half_length;
        if (extend_cnv_regions)
        {
            //printMessage("Querying preceding SNPs...");
            sv_half_length = (end_pos - start_pos) / 2;
            int64_t before_sv_start = start_pos - sv_half_length;
            std::pair<SNPData, bool> preceding_snps = this->querySNPRegion(chr, before_sv_start, start_pos-1, snp_info, pos_depth_map, mean_chr_cov);
            //printMessage("Running Viterbi algorithm for preceding SNPs...");
            std::vector<int> preceding_states = runViterbi(hmm, preceding_snps.first);
            //printMessage("Updating SNP data for preceding SNPs...");
            SNPData& pre_snp = preceding_snps.first;
            this->updateSNPVectors(snp_data, pre_snp.pos, pre_snp.pfb, pre_snp.baf, pre_snp.log2_cov, preceding_states, pre_snp.is_snp);
        }

        // Within SV SNPs:
        //printMessage("Updating SNP data for SV SNPs...");
        this->updateSNPVectors(snp_data, sv_snps.pos, sv_snps.pfb, sv_snps.baf, sv_snps.log2_cov, state_sequence, sv_snps.is_snp);

        // Following SNPs:
        if (extend_cnv_regions)
        {
            //printMessage("Querying following SNPs...");
            int64_t after_sv_end = end_pos + sv_half_length;
            std::pair<SNPData, bool> following_snps = this->querySNPRegion(chr, end_pos+1, after_sv_end, snp_info, pos_depth_map, mean_chr_cov);
            //printMessage("Running Viterbi algorithm for following SNPs...");
            std::vector<int> following_states = runViterbi(hmm, following_snps.first);
            //printMessage("Updating SNP data for following SNPs...");
            SNPData& post_snp = following_snps.first;
            this->updateSNPVectors(snp_data, post_snp.pos, post_snp.pfb, post_snp.baf, post_snp.log2_cov, following_states, post_snp.is_snp);
            //printMessage("Finished processing surrounding region for SV " + std::to_string(current_sv) + " of " + std::to_string(sv_count) + "...");
        }
    }

    // Print the SV type counts (proportions)
    printMessage("SV type counts for chromosome " + chr + ":");
    for (auto const& type_count : cnv_type_counts)
    {
        int type = type_count.first;
        int count = type_count.second;
        double pct = (double) count / (double) sv_count * 100;
        printMessage(sv_types::SVTypeString[type] + ": " + std::to_string(count) + " / " + std::to_string(sv_count) + " (" + std::to_string(pct) + "%)");
    }

    return snp_data;
}

void CNVCaller::updateSVType(std::map<SVCandidate, SVInfo> &sv_candidates, SVCandidate key, int sv_type, std::string data_type)
{
    // Lock the SV candidate map
    std::lock_guard<std::mutex> lock(this->sv_candidates_mtx);

    // Update the SV type if it is not unknown
    if (sv_type != sv_types::UNKNOWN)
    {
        sv_candidates[key].sv_type = sv_type;
        sv_candidates[key].data_type.insert(data_type);
    }
}

void CNVCaller::updateSVGenotype(std::map<SVCandidate,SVInfo>& sv_candidates, SVCandidate key, std::string genotype)
{
    // Lock the SV candidate map
    std::lock_guard<std::mutex> lock(this->sv_candidates_mtx);

    // Update the SV genotype
    sv_candidates[key].genotype = genotype;
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

    // Get the window size for HMM observations
    int window_size = this->input_data->getWindowSize();

    // Loop over SV call chromosomes
    std::cout << "Predicting copy number states for SV candidates..." << std::endl;
    std::set<std::string> sv_chromosomes = sv_calls.getChromosomes();
    // SNPData snp_data;
    int chr_count = (int) sv_chromosomes.size();
    int current_chr = 0;

    // Get the number of threads
    int num_threads = this->input_data->getThreadCount();

    // Process each chromosome in parallel
    int thread_count = 0;
    std::vector<std::future<SNPData>> futures;
    SNPData snp_data;
    for (auto const& chr : sv_chromosomes)
    {
        // Get the SNP information for the chromosome
        //std::cout << "Getting SNP information for chromosome " << chr << "..." << std::endl;

        // Get the SV candidates for the chromosome
        std::map<SVCandidate, SVInfo>& chr_sv_calls = sv_calls.getChromosomeSVs(chr);

        // Start a thread to run the copy number prediction for the entire chromosome
        //std::cout << "Running copy number prediction for chromosome " << chr
        //<< "..." << std::endl;
        printMessage("Running copy number prediction for chromosome " + chr + "...");

        // runCopyNumberPrediction(chr, sv_calls, snp_info, snp_data, hmm,
        // window_size);
        // Use futures to run the copy number prediction for the entire
        // chromosome and get the SNP data later
        std::future<SNPData> future = std::async(std::launch::async, &CNVCaller::runCopyNumberPrediction, this, chr, std::ref(chr_sv_calls), std::ref(snp_info), hmm, window_size);
        futures.push_back(std::move(future));
        thread_count++;

        // If the number of threads is reached, then wait for the threads to
        // finish and process the SNP data
        if (thread_count == num_threads)
        {
            // Loop through the futures and get the SNP data
            printMessage("Waiting for threads to finish...");
            for (auto& future : futures)
            {
                // Wait for the future to finish
                future.wait();

                // Get the SNP data from the future
                SNPData chr_snp_data = future.get();

                // Add the SNP data to the SNP data map
                printMessage("Updating SNP data...");
                // Lock the SNP data mutex
                this->snp_mtx.lock();
                updateSNPVectors(snp_data, chr_snp_data.pos, chr_snp_data.pfb, chr_snp_data.baf, chr_snp_data.log2_cov, chr_snp_data.state_sequence, chr_snp_data.is_snp);
                this->snp_mtx.unlock();

                // Update count and print progress
                current_chr++;
                std::cout << "Finished predicting copy number states for " << current_chr << " of " << chr_count << " chromosomes" << std::endl;
            }

            // Reset the thread count
            thread_count = 0;
            futures.clear();
        }
    }

    // Wait for the remaining threads to finish
    int remaining_threads = (int) futures.size();
    printMessage("Waiting for remaining " + std::to_string(remaining_threads) + " threads to finish...");
    for (auto& future : futures)
    {
        // Wait for the future to finish
        future.wait();

        // Get the SNP data from the future
        SNPData chr_snp_data = future.get();

        // Add the SNP data to the SNP data map
        printMessage("Updating SNP data...");
        // Lock the SNP data mutex
        this->snp_mtx.lock();
        updateSNPVectors(snp_data, chr_snp_data.pos, chr_snp_data.pfb, chr_snp_data.baf, chr_snp_data.log2_cov, chr_snp_data.state_sequence, chr_snp_data.is_snp);
        this->snp_mtx.unlock();

        // Update count and print progress
        current_chr++;
        std::cout << "Finished predicting copy number states for " << current_chr << " of " << chr_count << " chromosomes" << std::endl;
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

double CNVCaller::calculateLog2Ratio(int start_pos, int end_pos, std::unordered_map<uint64_t, int> &pos_depth_map, double mean_chr_cov)
{
    // // Print the size in kb if greater than 10 kb
    // int region_size = end_pos - start_pos;
    // int size_threshold = 10000;
    // if (region_size > size_threshold)
    // {
    //     std::cout << "Calculating log2 ratio, length (kb): " << region_size / 1000 << std::endl;
    // }

    // Use the position and depth map to calculate the log2 ratio
    double cum_depth = 0;
    int pos_count = 0;
    for (int i = start_pos; i <= end_pos; i++)
    {
        // Check if the position is in the map
        auto it = pos_depth_map.find(i);
        if (it == pos_depth_map.end())
        {
            //std::cout << "Position " << i << " not found in depth map" << std::endl;
            continue;
        }

        // Get the depth for the position
        int depth = pos_depth_map[i];

        // Update the position count and cumulative depth
        pos_count++;
        cum_depth += depth;
    }

    // Continue if there are no positions in the region
    if (pos_count == 0)
    {
        //std::cout << "No positions found in region " << start_pos << "-" << end_pos << std::endl;
        return 0;
    }

    // Calculate the window coverage log2 ratio
    double window_mean_cov = (double) cum_depth / (double) pos_count;

    // Calculate the log2 ratio for the window
    double window_log2_ratio = log2(window_mean_cov / mean_chr_cov);

    // if (region_size > size_threshold)
    // {
    //     std::cout << "Finished calculating large log2 ratio, length (kb): " << region_size / 1000 << std::endl;
    // }

    return window_log2_ratio;
}

void CNVCaller::readSNPAlleleFrequencies(std::string snp_filepath, SNPInfo &snp_info, bool whole_genome)
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
                // chromosome within the SNP region, filtering for SNPS only,
                // and within the MIN-MAX range of frequencies.
                std::string filter_criteria = "INFO/variant_type=\"snv\" && AF >= " + std::to_string(MIN_PFB) + " && AF <= " + std::to_string(MAX_PFB);
                std::string cmd = \
                    "bcftools query -r " + region_chunk + " -f '%POS\t%AF\n' -i '" + filter_criteria + "' " + pfb_filepath + " 2>/dev/null";
                    //"bcftools query -r " + region_chunk + " -f '%POS\t%AF\n' -i 'INFO/variant_type=\"snv\"' " + pfb_filepath + " 2>/dev/null";

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

void CNVCaller::updateSNPData(SNPData& snp_data, int64_t pos, double pfb, double baf, double log2_cov, bool is_snp)
{
    // Update the SNP data
    snp_data.pos.emplace_back(pos);
    snp_data.pfb.emplace_back(pfb);
    snp_data.baf.emplace_back(baf);
    snp_data.log2_cov.emplace_back(log2_cov);
    snp_data.is_snp.emplace_back(is_snp);
}

void CNVCaller::updateSNPVectors(SNPData &snp_data, std::vector<int64_t> &pos, std::vector<double> &pfb, std::vector<double> &baf, std::vector<double> &log2_cov, std::vector<int> &state_sequence, std::vector<bool> &is_snp)
{
    // Update the SNP data by first reserving space for the new data and then
    // inserting the new data
    snp_data.pos.reserve(snp_data.pos.size() + pos.size());
    snp_data.pos.insert(snp_data.pos.end(), std::make_move_iterator(pos.begin()), std::make_move_iterator(pos.end()));

    snp_data.pfb.reserve(snp_data.pfb.size() + pfb.size());
    snp_data.pfb.insert(snp_data.pfb.end(), std::make_move_iterator(pfb.begin()), std::make_move_iterator(pfb.end()));

    snp_data.baf.reserve(snp_data.baf.size() + baf.size());
    snp_data.baf.insert(snp_data.baf.end(), std::make_move_iterator(baf.begin()), std::make_move_iterator(baf.end()));

    snp_data.log2_cov.reserve(snp_data.log2_cov.size() + log2_cov.size());
    snp_data.log2_cov.insert(snp_data.log2_cov.end(), std::make_move_iterator(log2_cov.begin()), std::make_move_iterator(log2_cov.end()));

    snp_data.state_sequence.reserve(snp_data.state_sequence.size() + state_sequence.size());
    snp_data.state_sequence.insert(snp_data.state_sequence.end(), std::make_move_iterator(state_sequence.begin()), std::make_move_iterator(state_sequence.end()));

    snp_data.is_snp.reserve(snp_data.is_snp.size() + is_snp.size());
    snp_data.is_snp.insert(snp_data.is_snp.end(), std::make_move_iterator(is_snp.begin()), std::make_move_iterator(is_snp.end()));
}
