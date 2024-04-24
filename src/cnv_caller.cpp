
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
// std::max
#include <algorithm>


#include "utils.h"
#include "sv_data.h"
#include "sv_types.h"

#define MIN_PFB 0.01
#define MAX_PFB 0.99
/// @endcond

using namespace sv_types;

// Function to call the Viterbi algorithm for the CHMM
std::vector<int> CNVCaller::runViterbi(CHMM hmm, SNPData& snp_data)
{
    int data_count = (int) snp_data.pos.size();
    std::lock_guard<std::mutex> lock(this->hmm_mtx);  // Lock the mutex for the HMM
    std::vector<int> state_sequence = testVit_CHMM(hmm, data_count, snp_data.log2_cov, snp_data.baf, snp_data.pfb);
    return state_sequence;
}

// Function to obtain SNP information for a region
std::pair<SNPData, bool> CNVCaller::querySNPRegion(std::string chr, int64_t start_pos, int64_t end_pos, SNPInfo& snp_info, std::unordered_map<uint64_t, int>& pos_depth_map, double mean_chr_cov)
{
    // Get the SNP information for the region
    SNPData snp_data;
    bool snps_found = false;
    int window_size = this->input_data->getWindowSize();

    for (int64_t i = start_pos; i <= end_pos; i += window_size)
    {
        // Run a sliding non-overlapping window of size window_size across
        // the SV region and calculate the log2 ratio for each window
        int64_t window_start = i;
        int64_t window_end = std::min(i + window_size - 1, end_pos);

        // Get the SNP info for the window
        this->snp_data_mtx.lock();
        std::tuple<std::vector<int64_t>, std::vector<double>, std::vector<double>> window_snps = snp_info.querySNPs(chr, window_start, window_end);
        this->snp_data_mtx.unlock();
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
            double window_log2_ratio = calculateLog2Ratio(window_start, window_end, pos_depth_map, mean_chr_cov);

            // Use the window center position
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
                double bin_cov = calculateLog2Ratio(bin_start, bin_end, pos_depth_map, mean_chr_cov);
                this->updateSNPData(snp_data, (bin_start + bin_end) / 2, snp_window_pfbs[j], snp_window_bafs[j], bin_cov, true);

                // Update the previous bin start
                bin_start = bin_end + 1;
            }
        }
    }

    return std::make_pair(snp_data, snps_found);
}
SNPData CNVCaller::runCopyNumberPrediction(std::string chr, std::map<SVCandidate, SVInfo>& sv_candidates, SNPInfo& snp_info, CHMM hmm, int window_size, double mean_chr_cov)
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

    // If extending the CNV regions, then extend the SV region by window size *
    // N. Otherwise, log2 ratios will be zero due to missing read depth data
    // before/after the first/last SV positions
    if (this->input_data->getSaveCNVData())
    {
        int extend_factor = 100;
        first_pos = std::max((int64_t) 1, first_pos - window_size * extend_factor);
        last_pos = last_pos + window_size * extend_factor;
    }

    // Generate the read depth map for the entire SV region
    std::unordered_map<uint64_t, int> pos_depth_map;
    printMessage("Calculating read depths for SV region " + chr + ":" + std::to_string((int)first_pos) + "-" + std::to_string((int)last_pos) + "...");
    calculateDepthsForSNPRegion(chr, first_pos, last_pos, pos_depth_map);
    
    // Run copy number prediction for the SV candidates
    // Loop through each SV candidate and predict the copy number state
    printMessage("Predicting copy number states for chromosome " + chr + "...");

    // Create a map with counts for each CNV type
    std::map<int, int> cnv_type_counts;
    for (int i = 0; i < 6; i++)
    {
        cnv_type_counts[i] = 0;
    }

    // Split the SV candidates into chunks for each thread
    int chunk_count = this->input_data->getThreadCount();
    std::vector<std::vector<SVCandidate>> sv_chunks = splitSVCandidatesIntoChunks(sv_candidates, chunk_count);

    // Loop through each SV chunk and run the copy number prediction in parallel
    std::vector<std::future<SNPData>> futures;
    for (const auto& sv_chunk : sv_chunks)
    {
        // Run the copy number prediction for the SV chunk
        std::future<SNPData> future = std::async(std::launch::async, &CNVCaller::runCopyNumberPredictionChunk, this, chr, std::ref(sv_candidates), sv_chunk, std::ref(snp_info), hmm, window_size, mean_chr_cov, std::ref(pos_depth_map));
        futures.push_back(std::move(future));
    }

    // Get the SNP data for each SV chunk
    int current_chunk = 0;
    for (auto& future : futures)
    {
        current_chunk++;
        SNPData chunk_snp_data = std::move(future.get());
        if (this->input_data->getVerbose())
        {
            printMessage("Finished processing SV chunk " + std::to_string(current_chunk) + " of " + std::to_string(chunk_count) + "...");
        }

        // Update the SNP data
        if (this->input_data->getSaveCNVData())
        {
            this->updateSNPVectors(snp_data, chunk_snp_data.pos, chunk_snp_data.pfb, chunk_snp_data.baf, chunk_snp_data.log2_cov, chunk_snp_data.state_sequence, chunk_snp_data.is_snp);
            if (this->input_data->getVerbose())
            {
                printMessage("Updated SNP data for SV chunk " + std::to_string(current_chunk) + " of " + std::to_string(chunk_count) + "...");
            }
        }
    }

    printMessage("Finished predicting copy number states for chromosome " + chr + "...");

    return snp_data;
}

SNPData CNVCaller::runCopyNumberPredictionChunk(std::string chr, std::map<SVCandidate, SVInfo>& sv_candidates, std::vector<SVCandidate> sv_chunk, SNPInfo& snp_info, CHMM hmm, int window_size, double mean_chr_cov, std::unordered_map<uint64_t, int>& pos_depth_map)
{
    SNPData snp_data;

    // Map with counts for each CNV type
    std::map<int, int> cnv_type_counts;
    for (int i = 0; i < 6; i++)
    {
        cnv_type_counts[i] = 0;
    }
    
    // Loop through each SV candidate and predict the copy number state
    for (const auto& sv_call : sv_chunk)
    {
        // Get the SV candidate
        const SVCandidate& candidate = sv_call;

        // Get the start and end positions of the SV call
        int64_t start_pos = std::get<0>(candidate);
        int64_t end_pos = std::get<1>(candidate);

        // Get the depth at the start position. This is used as the FORMAT/DP
        // value in the VCF file
        int dp_value = pos_depth_map[start_pos];
        this->updateDPValue(sv_candidates, sv_call, dp_value);

        // Loop through the SV region, calculate the log2 ratios, and run the
        // Viterbi algorithm to predict the copy number states
        std::pair<SNPData, bool> snp_call = this->querySNPRegion(chr, start_pos, end_pos, snp_info, pos_depth_map, mean_chr_cov);
        SNPData& sv_snps = snp_call.first;
        bool snps_found = snp_call.second;

        // Run the Viterbi algorithm
        std::vector<int> state_sequence = runViterbi(hmm, sv_snps);

        // Determine if there is a majority state and if it is greater than 75%
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
        double pct_threshold = 0.75;
        int state_count = (int) state_sequence.size();
        if ((double) max_count / (double) state_count < pct_threshold)
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

        // Update the SV type if it is not unknown
        if (cnv_type != sv_types::UNKNOWN)
        {
            this->updateSVType(sv_candidates, sv_call, cnv_type, data_type);

            // Update the SV type counts
            cnv_type_counts[cnv_type]++;
        }

        // Update the SV genotype
        this->updateSVGenotype(sv_candidates, sv_call, genotype);

        if (this->input_data->getSaveCNVData())
        {
            // Preceding SNPs:
            int64_t sv_half_length = (end_pos - start_pos) / 2;
            int64_t before_sv_start = start_pos - sv_half_length;
            std::pair<SNPData, bool> preceding_snps = this->querySNPRegion(chr, before_sv_start, start_pos-1, snp_info, pos_depth_map, mean_chr_cov);
            std::vector<int> preceding_states = runViterbi(hmm, preceding_snps.first);
            SNPData& pre_snp = preceding_snps.first;
            this->updateSNPVectors(snp_data, pre_snp.pos, pre_snp.pfb, pre_snp.baf, pre_snp.log2_cov, preceding_states, pre_snp.is_snp);

            // Within SV SNPs:
            this->updateSNPVectors(snp_data, sv_snps.pos, sv_snps.pfb, sv_snps.baf, sv_snps.log2_cov, state_sequence, sv_snps.is_snp);

            // Following SNPs:
            int64_t after_sv_end = end_pos + sv_half_length;
            std::pair<SNPData, bool> following_snps = this->querySNPRegion(chr, end_pos+1, after_sv_end, snp_info, pos_depth_map, mean_chr_cov);
            std::vector<int> following_states = runViterbi(hmm, following_snps.first);
            SNPData& post_snp = following_snps.first;
            this->updateSNPVectors(snp_data, post_snp.pos, post_snp.pfb, post_snp.baf, post_snp.log2_cov, following_states, post_snp.is_snp);
        }
    }

    // Print the SV type counts (proportions)
    if (this->input_data->getVerbose())
    {
        printMessage("SV chunk type counts for chromosome " + chr + ":");
        for (auto const& type_count : cnv_type_counts)
        {
            int type = type_count.first;
            int count = type_count.second;
            double pct = (double) count / (double) sv_chunk.size() * 100;
            printMessage(sv_types::SVTypeString[type] + ": " + std::to_string(count) + " / " + std::to_string(sv_chunk.size()) + " (" + std::to_string(pct) + "%)");
        }
    }

    return snp_data;
}

void CNVCaller::updateSVType(std::map<SVCandidate, SVInfo> &sv_candidates, SVCandidate key, int sv_type, std::string data_type)
{
    // Lock the SV candidate map
    std::lock_guard<std::mutex> lock(this->sv_candidates_mtx);

    // Update the SV type if the update is not unknown
    if (sv_type != sv_types::UNKNOWN)
    {
        // Update the SV type if the existing type is unknown
        if (sv_candidates[key].sv_type == sv_types::UNKNOWN)
        {
            sv_candidates[key].sv_type = sv_type;
        }

        // sv_candidates[key].sv_type = sv_type;
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

void CNVCaller::updateDPValue(std::map<SVCandidate,SVInfo>& sv_candidates, SVCandidate key, int dp_value)
{
    // Lock the SV candidate map
    std::lock_guard<std::mutex> lock(this->sv_candidates_mtx);

    // Update the DP value
    sv_candidates[key].read_depth = dp_value;
}

std::vector<std::string> CNVCaller::splitRegionIntoChunks(std::string chr, int64_t start_pos, int64_t end_pos, int chunk_count)
{
    // Split the region into chunks
    std::vector<std::string> region_chunks;
    int64_t region_length = end_pos - start_pos + 1;
    int64_t chunk_size = std::ceil((double) region_length / (double) chunk_count);
    int64_t chunk_start = start_pos;
    int64_t chunk_end = 0;
    for (int i = 0; i < chunk_count; i++)
    {
        chunk_end = chunk_start + chunk_size - 1;

        if (i == chunk_count - 1)
        {
            chunk_end = end_pos;
        }

        // Add the region chunk to the vector
        region_chunks.push_back(chr + ":" + std::to_string(chunk_start) + "-" + std::to_string(chunk_end));

        // Update the chunk start
        chunk_start = chunk_end + 1;
    }

    return region_chunks;
}

std::vector<std::vector<SVCandidate>> CNVCaller::splitSVCandidatesIntoChunks(std::map<SVCandidate,SVInfo>& sv_candidates, int chunk_count)
{
    // Split the SV candidates into chunks
    std::vector<std::vector<SVCandidate>> sv_chunks;
    int sv_count = (int) sv_candidates.size();
    int chunk_size = std::ceil((double) sv_count / (double) chunk_count);
    int current_chunk = 0;
    std::vector<SVCandidate> current_sv_chunk;
    for (auto const& sv_call : sv_candidates)
    {
        // Add the SV candidate to the current chunk
        current_sv_chunk.push_back(sv_call.first);

        // If the current chunk size is reached, then add the chunk to the
        // vector and reset the current chunk
        if ((int) current_sv_chunk.size() == chunk_size)
        {
            sv_chunks.push_back(current_sv_chunk);
            current_sv_chunk.clear();
            current_chunk++;
        }
    }

    // Add the remaining SV candidates to the last chunk
    if (current_sv_chunk.size() > 0)
    {
        sv_chunks.push_back(current_sv_chunk);
    }

    return sv_chunks;
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

    // Read the HMM from file
    std::string hmm_filepath = this->input_data->getHMMFilepath();
    std::cout << "Reading HMM from file: " << hmm_filepath << std::endl;
    CHMM hmm = ReadCHMM(hmm_filepath.c_str());

    // Get the window size for HMM observations
    int window_size = this->input_data->getWindowSize();

    // Loop over SV call chromosomes
    std::cout << "Predicting copy number states for SV candidates..." << std::endl;
    std::set<std::string> sv_chromosomes = sv_calls.getChromosomes();

    // Process each chromosome
    SNPData snp_data;
    int current_chr = 0;
    int chr_count = (int) sv_chromosomes.size();
    for (auto const& chr : sv_chromosomes)
    {
        // Get the SV candidates for the chromosome
        std::map<SVCandidate, SVInfo>& chr_sv_calls = sv_calls.getChromosomeSVs(chr);

        // First, calculate the mean chromosome coverage
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

        // Read the SNP positions and B-allele frequency values from the VCF file
        SNPInfo snp_info;
        std::cout << "Reading SNP allele frequencies for chromosome " << chr << " from VCF file..." << std::endl;
        std::string snp_filepath = this->input_data->getSNPFilepath();
        readSNPAlleleFrequencies(chr, snp_filepath, snp_info);

        // Get the population frequencies for each SNP
        std::cout << "Obtaining SNP population frequencies for chromosome " << chr << "..." << std::endl;
        getSNPPopulationFrequencies(chr, snp_info);

        // Start a thread to run the copy number prediction for the entire chromosome
        printMessage("Running copy number prediction for chromosome " + chr + "...");

        // Run copy number prediction for the chromosome
        SNPData chr_snp_data = runCopyNumberPrediction(chr, chr_sv_calls, snp_info, hmm, window_size, mean_chr_cov);

        // Add the SNP data to the SNP data map
        if (this->input_data->getSaveCNVData())
        {
            printMessage("Updating SNP data...");
            updateSNPVectors(snp_data, chr_snp_data.pos, chr_snp_data.pfb, chr_snp_data.baf, chr_snp_data.log2_cov, chr_snp_data.state_sequence, chr_snp_data.is_snp);
        }
        printMessage("Finished processing chromosome " + chr + " (" + std::to_string(++current_chr) + " of " + std::to_string(chr_count) + ")...");
    }

    std::cout << "Finished predicting copy number states for all chromosomes" << std::endl;
    std::cout << "Found a total of " << sv_calls.totalCalls() << " SV candidates" << std::endl;

    // Save the SNP data to a TSV file
    if (this->input_data->getSaveCNVData())
    {
        std::string snp_output_tsv = this->input_data->getOutputDir() + "/cnv_data.tsv";
        saveToTSV(snp_data, snp_output_tsv);
        std::cout << "Saved SNP data to " << snp_output_tsv << std::endl;
    }
}

// Calculate the mean chromosome coverage
double CNVCaller::calculateMeanChromosomeCoverage(std::string chr)
{
    std::string input_filepath = this->input_data->getShortReadBam();

    // Get the number of threads
    int num_threads = this->input_data->getThreadCount();

    // Split the chromosome into equal parts for each thread
    int chr_len = this->input_data->getRefGenomeChromosomeLength(chr);

    // Split the chromosome into equal parts for each thread
    std::vector<std::string> region_chunks = splitRegionIntoChunks(chr, 1, chr_len, num_threads);

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
                printError("ERROR: Could not open pipe for command: " + std::string(cmd));
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
        std::tuple<uint64_t, uint64_t> result = std::move(future.get());

        // Update the position count and cumulative depth
        pos_count += std::get<0>(result);
        cum_depth += std::get<1>(result);
    }

    // Calculate the mean chromosome coverage
    double mean_chr_cov = (double) cum_depth / (double) pos_count;

    return mean_chr_cov;
}

void CNVCaller::calculateDepthsForSNPRegion(std::string chr, int64_t start_pos, int64_t end_pos, std::unordered_map<uint64_t, int>& pos_depth_map)
{
    std::string input_filepath = this->input_data->getShortReadBam();

    // Get the number of threads
    int num_threads = this->input_data->getThreadCount();

    // Split the region into equal parts for each thread if the region is larger
    // than 10 kb
    std::vector<std::string> region_chunks;
    int64_t region_size = end_pos - start_pos;
    int64_t min_threading_size = 10000;
    printMessage("Region size: " + std::to_string(region_size));
    if (region_size < min_threading_size)
    {
        region_chunks.push_back(chr + ":" + std::to_string(start_pos) + "-" + std::to_string(end_pos));
    } else {
        region_chunks = splitRegionIntoChunks(chr, start_pos, end_pos, num_threads);
    }

    // Loop through each region chunk and get the mean chromosome coverage in
    // parallel
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
    int current_chunk = 0;
    for (auto& future : futures)
    {
        current_chunk++;

        // Wait for the future to finish
        future.wait();

        // Get the result from the future
        std::unordered_map<uint64_t, int> result = std::move(future.get());

        // Merge the result into the map of positions and depths
        this->mergePosDepthMaps(pos_depth_map, result);
        if (this->input_data->getVerbose())
        {
            printMessage("Completed region chunk " + std::to_string(current_chunk) + " of " + std::to_string(region_chunks.size()) + "...");
        }
    }
}

void CNVCaller::mergePosDepthMaps(std::unordered_map<uint64_t, int>& pos_depth_map, std::unordered_map<uint64_t, int>& result)
{ 
    // Reserve space for the new map
    pos_depth_map.reserve(pos_depth_map.size() + result.size());

    // Move elements from the new map to the existing map
    for (auto& elem : result)
    {
        // Use the move assignment operator to move the element from the new
        // map to the existing map
        pos_depth_map[elem.first] = std::move(elem.second);
    }
}

double CNVCaller::calculateLog2Ratio(int start_pos, int end_pos, std::unordered_map<uint64_t, int> &pos_depth_map, double mean_chr_cov)
{
    // Use the position and depth map to calculate the log2 ratio
    double cum_depth = 0;
    int pos_count = 0;
    for (int i = start_pos; i <= end_pos; i++)
    {
        // Check if the position is in the map
        auto it = pos_depth_map.find(i);
        if (it == pos_depth_map.end())
        {
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
        return 0;
    }

    // Calculate the window coverage log2 ratio
    double window_mean_cov = (double) cum_depth / (double) pos_count;

    // Calculate the log2 ratio for the window
    double window_log2_ratio = log2(window_mean_cov / mean_chr_cov);

    return window_log2_ratio;
}

void CNVCaller::readSNPAlleleFrequencies(std::string chr, std::string filepath, SNPInfo& snp_info)
{
    // Create a VCF filepath of filtered SNPs
    std::string filtered_snp_vcf_filepath = this->input_data->getOutputDir() + "/filtered_snps.vcf";

    // Check that the SNP file is sorted by running bcftools index and reading
    // the error output
    std::string index_cmd = "bcftools index " + filepath + " 2>&1 | grep -i error";
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

    // Filter variants by depth, quality, and region
    if (this->input_data->getVerbose()) {
        std::cout << "Filtering SNPs by depth, quality, and region..." << std::endl;
    }
    std::string cmd = "bcftools view -r " + chr + " -v snps -i 'QUAL > 30 && DP > 10 && FILTER = \"PASS\"' " + filepath + " > " + filtered_snp_vcf_filepath;
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
    if (this->input_data->getVerbose()) {
        std::cout << "Extracting B-allele frequency data from filtered SNPs..." << std::endl;
    }
    cmd = "bcftools query -f '%POS,[%AD]\n' " + filtered_snp_vcf_filepath + " 2>/dev/null";
    FILE *fp = popen(cmd.c_str(), "r");
    if (fp == NULL)
    {
        std::cerr << "ERROR: Could not open pipe for command: " << cmd << std::endl;
        exit(1);
    }

    // Read the reference and alternate allele depths from the VCF file
    std::string alt_allele = "";  // Alternate allele
    uint64_t pos = 0;
    int ref_ad = 0;
    int alt_ad = 0;
    const int line_size = 256;
    char line[line_size];  // Line buffer
    std::vector<int64_t> locations;
    std::vector<double> bafs;
    while (fgets(line, line_size, fp) != NULL)
    {
        // Parse the line
        char *tok = strtok(line, ",");  // Tokenize the line
        int col = 0;  // Column index
        while (tok != NULL)
        {
            // Get the position from column 2
            if (col == 0)
            {
                pos = atoi(tok);
            }

            // Get the AD for the reference allele from column 3
            else if (col == 1)
            {
                ref_ad = atoi(tok);
            }

            // Get the AD for the non-reference allele from column 4
            else if (col == 2)
            {
                alt_ad = atoi(tok);
            }

            // Move to the next token
            tok = strtok(NULL, ",");
            col++;
        }

        // Calculate the B-allele frequency (BAF) as the ratio of the alternate
        // allele depth to the total depth (reference + alternate)
        double baf = (double) alt_ad / (double) (ref_ad + alt_ad);

        // Add a new location and BAF value to the chromosome's SNP data
        // (population frequency and log2 ratio will be added later)
        snp_info.insertSNPAlleleFrequency(chr, pos, baf);
    }

    // Close the pipe
    pclose(fp);

    if (this->input_data->getVerbose()) {
        std::cout << "Finished extracting B-allele frequency data from filtered SNPs" << std::endl;
    }
}

void CNVCaller::getSNPPopulationFrequencies(std::string chr, SNPInfo& snp_info)
{
    // Get the population frequency file for the chromosome
    std::string pfb_filepath = this->input_data->getAlleleFreqFilepath(chr);

    // Determine the ethnicity-specific allele frequency key
    std::string AF_key = "AF";
    if (this->input_data->getEthnicity() != "")
    {
        AF_key += "_" + this->input_data->getEthnicity();
    }

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
        return;
    }

    std::cout << "Reading population frequencies for chromosome " << chr << " from " << pfb_filepath << std::endl;

    // Get the start and end SNP positions for the chromosome (1-based
    // index)
    std::pair<int64_t, int64_t> snp_range = snp_info.getSNPRange(chr);
    int64_t snp_start = snp_range.first;
    int64_t snp_end = snp_range.second;
    std::cout << "SNP range for chromosome " << chr << ": " << snp_start << "-" << snp_end << std::endl;

    // Get the number of avaiable threads
    int num_threads = this->input_data->getThreadCount();

    // Get the region chunks
    std::vector<std::string> region_chunks = splitRegionIntoChunks(chr_gnomad, snp_start, snp_end, num_threads);

    // Loop through each region chunk and get the population frequencies in
    // parallel
    std::unordered_map<int, double> pos_pfb_map;
    std::vector<std::thread> threads;
    std::vector<std::future<std::unordered_map<int, double>>> futures;
    for (const auto& region_chunk : region_chunks)
    {
        // Create a lambda function to get the population frequencies for the
        // region chunk
        auto get_pfb = [region_chunk, pfb_filepath, AF_key]() -> std::unordered_map<int, double>
        {
            // Run bcftools query to get the population frequencies for the
            // chromosome within the SNP region, filtering for SNPS only,
            // and within the MIN-MAX range of frequencies.
            // TODO: Update to use ethnicity-specific population frequencies
            // Example from gnomAD:
            // ##INFO=<ID=AF_asj,Number=A,Type=Float,Description="Alternate
            // allele frequency in samples of Ashkenazi Jewish ancestry">
            // std::string ethnicity_suffix = "_asj";  // Ashkenazi Jewish
            // (leave empty for all populations)
            std::string filter_criteria = "INFO/variant_type=\"snv\" && " + AF_key + " >= " + std::to_string(MIN_PFB) + " && " + AF_key + " <= " + std::to_string(MAX_PFB);
            std::string cmd = \
                "bcftools query -r " + region_chunk + " -f '%POS\t%" + AF_key + "\n' -i '" + filter_criteria + "' " + pfb_filepath + " 2>/dev/null";

            // [TEST] Print the command
            std::cout << "Command: " << cmd << std::endl;

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
    int pfb_count = 0;
    for (auto& future : futures)
    {
        // Wait for the future to finish
        future.wait();

        // Get the result from the future
        std::unordered_map<int, double> result = std::move(future.get());

        // Loop through the result and add to SNPInfo
        // printMessage("Adding population frequencies to SNPInfo...");
        for (auto& pair : result)
        {
            int pos = pair.first;
            double pfb = pair.second;

            // Lock the SNPInfo mutex
            this->snp_data_mtx.lock();
            snp_info.insertSNPPopulationFrequency(chr_snp, pos, pfb);
            this->snp_data_mtx.unlock();

            // Increment the population frequency count
            pfb_count++;

            // [TEST] Print 15 values
            if (pfb_count < 15)
            {
                printMessage("Population frequency for " + chr + ":" + std::to_string(pos) + " = " + std::to_string(pfb));
            }
        }
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

