
#include "cnv_caller.h"

#include <htslib/sam.h>

/// @cond
#include <iostream>
#include <sstream>
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
#include <algorithm>  // std::max

#include "utils.h"
#include "sv_data.h"
#include "sv_types.h"

#define MIN_PFB 0.01
#define MAX_PFB 0.99
/// @endcond

using namespace sv_types;

// Function to call the Viterbi algorithm for the CHMM
std::pair<std::vector<int>, double> CNVCaller::runViterbi(CHMM hmm, SNPData& snp_data)
{
    int data_count = (int) snp_data.pos.size();
    if (data_count == 0)
    {
        throw std::runtime_error("Error: No SNP data found for Viterbi algorithm.");
    }
    std::lock_guard<std::mutex> lock(this->hmm_mtx);  // Lock the mutex for the HMM
    std::pair<std::vector<int>, double> state_sequence = testVit_CHMM(hmm, data_count, snp_data.log2_cov, snp_data.baf, snp_data.pfb);
    return state_sequence;
}

// Function to obtain SNP information for a region
std::pair<SNPData, bool> CNVCaller::querySNPRegion(std::string chr, int64_t start_pos, int64_t end_pos, SNPInfo& snp_info, std::unordered_map<uint32_t, int>& pos_depth_map, double mean_chr_cov)
{
    SNPData snp_data;
    bool snps_found = false;
    int window_size = this->input_data->getWindowSize();

    // printMessage("Querying SNPs for region " + chr + ":" + std::to_string(start_pos) + "-" + std::to_string(end_pos) + "...");
    for (int64_t i = start_pos; i <= end_pos; i += window_size)
    {
        // Run a sliding non-overlapping window of size window_size across
        // the SV region and calculate the log2 ratio for each window
        int64_t window_start = i;
        int64_t window_end = std::min(i + window_size - 1, end_pos);

        // Get the SNP info for the window
        // std::cout << "Querying SNPs for window " << chr << ":" << window_start << "-" << window_end << "..." << std::endl;
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
            double window_log2_ratio = calculateLog2Ratio(window_start, window_end, pos_depth_map, mean_chr_cov);
            double pfb_default = 0.5;
            double baf_default = -1.0;  // Use -1.0 to indicate no BAF data
            this->updateSNPData(snp_data, (window_start + window_end) / 2, pfb_default, baf_default, window_log2_ratio, false);

        } else {
            snps_found = true;

            // Loop through the SNPs and calculate the log2 ratios
            int64_t bin_start = window_start;
            int64_t bin_end = 0;
            for (int j = 0; j < snp_count; j++)
            {
                // SNP bin starts at 1/2 the distance between the previous SNP
                // and the current SNP, and ends at 1/2 the distance between
                // the current SNP and the next SNP. For the first SNP, the
                // bin starts at the window start and ends at 1/2 the distance
                // between the first SNP and the next SNP, and for the last
                // SNP, the bin starts at 1/2 the distance between the previous
                // SNP and the last SNP and ends at the window end.
                int64_t snp_pos = snp_window_pos[j];
                bin_end = snp_pos + (j == snp_count-1 ? (window_end - snp_pos) / 2 : (snp_window_pos[j+1] - snp_pos) / 2);

                // Calculate the log2 ratio for the SNP bin
                double bin_cov = calculateLog2Ratio(bin_start, bin_end, pos_depth_map, mean_chr_cov);
                this->updateSNPData(snp_data, snp_pos, snp_window_pfbs[j], snp_window_bafs[j], bin_cov, true);

                // Update the previous bin start
                bin_start = bin_end + 1;
            }
        }
    }

    return std::make_pair(snp_data, snps_found);
}

std::tuple<double, SVType, std::string, bool> CNVCaller::runCopyNumberPrediction(std::string chr, SVCandidate& candidate)
{
    // std::cout << "Running copy number prediction for SV pair " << chr << ":" << std::get<0>(sv_one) << "-" << std::get<1>(sv_one) << " and " << std::get<0>(sv_two) << "-" << std::get<1>(sv_two) << "..." << std::endl;
     // Get the start and end positions of the SV call
    int64_t start_pos = std::get<0>(candidate);
    int64_t end_pos = std::get<1>(candidate);

    // Run the Viterbi algorithm on SNPs in the SV region +/- 1/2
    // the SV length
    int64_t sv_length = (end_pos - start_pos) / 2.0;
    int64_t snp_start_pos = std::max((int64_t) 1, start_pos - sv_length);
    int64_t snp_end_pos = end_pos + sv_length;

    // Query the SNP region for the SV candidate
    std::pair<SNPData, bool> snp_call = querySNPRegion(chr, snp_start_pos, snp_end_pos, this->snp_info, this->pos_depth_map, this->mean_chr_cov);
    SNPData& sv_snps = snp_call.first;
    bool sv_snps_found = snp_call.second;

    // Run the Viterbi algorithm
    std::pair<std::vector<int>, double> prediction = runViterbi(this->hmm, sv_snps);
    std::vector<int>& state_sequence = prediction.first;
    double likelihood = prediction.second;

    // Get all the states in the SV region
    std::vector<int> sv_states;
    for (size_t i = 0; i < state_sequence.size(); i++)
    {
        if (sv_snps.pos[i] >= start_pos && sv_snps.pos[i] <= end_pos)
        {
            sv_states.push_back(state_sequence[i]);
        }
    }

    // Determine if there is a majority state within the SV region and if it
    // is greater than 75%
    double pct_threshold = 0.75;
    int max_state = 0;
    int max_count = 0;
    for (int i = 0; i < 6; i++)
    {
        int state_count = std::count(sv_states.begin(), sv_states.end(), i+1);
        if (state_count > max_count)
        {
            max_state = i+1;
            max_count = state_count;
        }
    }
    
    // Update SV type and genotype based on the majority state
    SVType predicted_cnv_type = SVType::UNKNOWN;
    std::string genotype = "./.";
    int state_count = (int) sv_states.size();
    if ((double) max_count / (double) state_count > pct_threshold)
    {
        predicted_cnv_type = getSVTypeFromCNState(max_state);
        genotype = cnv_genotype_map[max_state];
    }
    sv_snps.state_sequence = std::move(state_sequence);  // Move the state sequence to the SNP data

    // Save the SV calls as a TSV file if enabled
    bool copy_number_change = (predicted_cnv_type != SVType::UNKNOWN && predicted_cnv_type != SVType::NEUTRAL);
    if (this->input_data->getSaveCNVData() && copy_number_change && (end_pos - start_pos) > 10000)
    {
        std::string cnv_type_str = getSVTypeString(predicted_cnv_type);
        std::string sv_filename = this->input_data->getOutputDir() + "/" + cnv_type_str + "_" + chr + "_" + std::to_string((int) start_pos) + "-" + std::to_string((int) end_pos) + "_SPLITALN.tsv";
        std::cout << "Saving SV split-alignment copy number predictions to " << sv_filename << std::endl;
        this->saveSVCopyNumberToTSV(sv_snps, sv_filename, chr, start_pos, end_pos, cnv_type_str, likelihood);
        // this->saveSVCopyNumberToTSV(best_snp_data, sv_filename, chr, best_pos.first, best_pos.second, cnv_type_str, best_likelihood);
    }

    return std::make_tuple(likelihood, predicted_cnv_type, genotype, sv_snps_found);
}


SNPData CNVCaller::runCIGARCopyNumberPrediction(std::string chr, std::map<SVCandidate, SVInfo> &sv_candidates, int min_length)
{
    SNPInfo& snp_info = this->snp_info;
    CHMM& hmm = this->hmm;
    int window_size = this->input_data->getWindowSize();
    double mean_chr_cov = this->mean_chr_cov;
    SNPData snp_data;

    // Filter the SV candidates by length
    std::map<SVCandidate, SVInfo> filtered_sv_candidates;
    for (const auto& sv_call : sv_candidates)
    {
        int64_t start_pos = std::get<0>(sv_call.first);
        int64_t end_pos = std::get<1>(sv_call.first);
        if ((end_pos - start_pos) >= min_length)
        {
            filtered_sv_candidates[sv_call.first] = sv_call.second;
        }
    }
    sv_candidates = std::move(filtered_sv_candidates);
    int sv_count = (int) sv_candidates.size();
    if (sv_count == 0)
    {
        return snp_data;
    }

   
    printMessage("Predicting CIGAR string copy number states for chromosome " + chr + "...");

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
        std::async(std::launch::async, &CNVCaller::runCIGARCopyNumberPredictionChunk, this, chr, std::ref(sv_candidates), sv_chunk, std::ref(snp_info), hmm, window_size, mean_chr_cov, std::ref(this->pos_depth_map));
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
    }

    printMessage("Finished predicting copy number states for chromosome " + chr + "...");

    return snp_data;
}

void CNVCaller::runCIGARCopyNumberPredictionChunk(std::string chr, std::map<SVCandidate, SVInfo>& sv_candidates, std::vector<SVCandidate> sv_chunk, SNPInfo& snp_info, CHMM hmm, int window_size, double mean_chr_cov, std::unordered_map<uint32_t, int>& pos_depth_map)
{
    // printMessage("Running copy number prediction for " + std::to_string(sv_chunk.size()) + " SV candidates on chromosome " + chr + "...");
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
        int64_t start_pos = std::get<0>(candidate);
        int64_t end_pos = std::get<1>(candidate);

        // Get the depth at the start position. This is used as the FORMAT/DP
        // value in the VCF file
        int dp_value = pos_depth_map[start_pos];
        this->updateDPValue(sv_candidates, sv_call, dp_value);

        // Loop through the SV region +/- 1/2 SV length and run copy number
        // predictions
        int64_t sv_half_length = (end_pos - start_pos) / 2.0;
        int64_t query_start = std::max((int64_t) 1, start_pos - sv_half_length);
        int64_t query_end = end_pos + sv_half_length;

        // printMessage("Querying SNPs for SV " + chr + ":" +
        // std::to_string((int)start_pos) + "-" + std::to_string((int)end_pos) +
        // "...");
        std::pair<SNPData, bool> snp_call = this->querySNPRegion(chr, query_start, query_end, snp_info, pos_depth_map, mean_chr_cov);
        SNPData& sv_snps = snp_call.first;
        bool snps_found = snp_call.second;

        // Run the Viterbi algorithm
        std::pair<std::vector<int>, double> prediction = runViterbi(hmm, sv_snps);
        std::vector<int>& state_sequence = prediction.first;
        double likelihood = prediction.second;

        // Get all the states in the SV region
        std::vector<int> sv_states;
        for (size_t i = 0; i < state_sequence.size(); i++)
        {
            if (sv_snps.pos[i] >= start_pos && sv_snps.pos[i] <= end_pos)
            {
                sv_states.push_back(state_sequence[i]);
            }
        }

        // Determine if there is a majority state within the SV region and if it
        // is greater than 75%
        int max_state = 0;
        int max_count = 0;
        for (int i = 0; i < 6; i++)
        {
            int state_count = std::count(sv_states.begin(), sv_states.end(), i+1);
            if (state_count > max_count)
            {
                max_state = i+1;
                max_count = state_count;
            }
        }

        // If there is no majority state, then set the state to unknown
        double pct_threshold = 0.75;
        int state_count = (int) sv_states.size();
        if ((double) max_count / (double) state_count < pct_threshold)
        {
            max_state = 0;
        }

        // Update the SV calls with the CNV type and genotype
        SVType cnv_type = getSVTypeFromCNState(max_state);
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

        // Update the SV copy number data
        this->updateSVCopyNumber(sv_candidates, sv_call, cnv_type, data_type, genotype, likelihood);

        // Save the SV calls as a TSV file if enabled, if the SV type is
        // known, and the length is greater than 10 kb
        SVType updated_sv_type = sv_candidates[sv_call].sv_type;
        if (this->input_data->getSaveCNVData() && updated_sv_type != SVType::UNKNOWN && (end_pos - start_pos) > 10000)
        {
            // Add the state sequence to the SNP data (avoid copying the data)
            sv_snps.state_sequence = std::move(state_sequence);

            // Save the SV calls as a TSV file
            std::string cnv_type_str = getSVTypeString(updated_sv_type);
            std::string sv_filename = this->input_data->getOutputDir() + "/" + cnv_type_str + "_" + chr + "_" + std::to_string((int) start_pos) + "-" + std::to_string((int) end_pos) + "_CIGAR.tsv";
            // std::cout << "Saving SV CIGAR copy number predictions to " <<
            // sv_filename << std::endl;
            printMessage("Saving SV CIGAR copy number predictions to " + sv_filename);
            this->saveSVCopyNumberToTSV(sv_snps, sv_filename, chr, start_pos, end_pos, cnv_type_str, likelihood);
        }
    }
}

void CNVCaller::updateSVCopyNumber(std::map<SVCandidate, SVInfo> &sv_candidates, SVCandidate key, SVType sv_type_update, std::string data_type, std::string genotype, double hmm_likelihood)
{
    // Update SV data from the HMM copy number prediction
    // Lock the SV candidate map
    std::lock_guard<std::mutex> lock(this->sv_candidates_mtx);

    // Update the SV type if the update is not unknown, and if the types don't
    // conflict (To avoid overwriting previous calls)
    SVType current_sv_type = sv_candidates[key].sv_type;
    if ((sv_type_update != SVType::UNKNOWN) && ((current_sv_type == sv_type_update) || (current_sv_type == SVType::UNKNOWN)))
    {
        sv_candidates[key].sv_type = sv_type_update;  // Update the SV type
        sv_candidates[key].data_type.insert(data_type);  // Update the data type

        // Update the likelihood if it is greater than the existing likelihood,
        // or if it is currently unknown (0.0)
        double previous_likelihood = sv_candidates[key].hmm_likelihood;
        if (previous_likelihood == 0.0 || hmm_likelihood > previous_likelihood)
        {
            sv_candidates[key].hmm_likelihood = hmm_likelihood;
        }

        // Update the genotype
        sv_candidates[key].genotype = genotype;
    }
}

void CNVCaller::updateDPValue(std::map<SVCandidate,SVInfo>& sv_candidates, SVCandidate key, int dp_value)
{
    std::lock_guard<std::mutex> lock(this->sv_candidates_mtx);
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

void CNVCaller::loadChromosomeData(std::string chr)
{
    std::string hmm_filepath = this->input_data->getHMMFilepath();
    std::cout << "Reading HMM from file: " << hmm_filepath << std::endl;
    this->hmm = ReadCHMM(hmm_filepath.c_str());

    printMessage("Calculating mean chromosome coverage for " + chr + "...");
    mean_chr_cov = calculateMeanChromosomeCoverage(chr);
    printMessage("Mean chromosome coverage for " + chr + ": " + std::to_string(mean_chr_cov));
    this->mean_chr_cov = mean_chr_cov;

    std::cout << "Reading SNP allele frequencies for chromosome " << chr << " from VCF file..." << std::endl;
    std::string snp_filepath = this->input_data->getSNPFilepath();
    readSNPAlleleFrequencies(chr, snp_filepath, this->snp_info);

    std::cout << "Obtaining SNP population frequencies for chromosome " << chr << "..." << std::endl;
    getSNPPopulationFrequencies(chr, this->snp_info);
    std::cout << "Finished loading chromosome data for " << chr << std::endl;
}

// Calculate the mean chromosome coverage
double CNVCaller::calculateMeanChromosomeCoverage(std::string chr)
{
    // Split the chromosome into equal parts for each thread
    int num_threads = this->input_data->getThreadCount();
    uint32_t chr_len = this->input_data->getRefGenomeChromosomeLength(chr);
    std::vector<std::string> region_chunks = splitRegionIntoChunks(chr, 1, chr_len, num_threads);

    // Calculate the mean chromosome coverage in parallel
    uint32_t pos_count = 0;
    uint64_t cum_depth = 0;
    std::vector<std::future<std::tuple<uint32_t, uint32_t, std::unordered_map<uint32_t, int>>>> futures;
    std::string input_filepath = this->input_data->getShortReadBam();
    for (const auto& region_chunk : region_chunks)
    {
        // Create a lambda function to get the mean chromosome coverage for the
        // region chunk
        auto get_mean_chr_cov = [region_chunk, input_filepath]() -> std::tuple<uint32_t, uint32_t, std::unordered_map<uint32_t, int>>
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
                printError("ERROR: Could not open pipe for command: " + std::string(cmd));
                exit(EXIT_FAILURE);
            }

            // Parse the outputs (position and depth)
            std::unordered_map<uint32_t, int> pos_depth_map;
            const int line_size = 256;
            char line[line_size];
            uint32_t pos;
            int depth;
            uint32_t pos_count = 0;
            uint64_t cum_depth = 0;
            while (fgets(line, line_size, fp) != NULL)
            {
                if (sscanf(line, "%u%d", &pos, &depth) == 2)
                {
                    pos_depth_map[pos] = depth;
                    pos_count++;
                    cum_depth += depth;
                }
            }
            pclose(fp);  // Close the process

            return std::make_tuple(pos_count, cum_depth, pos_depth_map);
        };
        std::future<std::tuple<uint32_t, uint32_t, std::unordered_map<uint32_t, int>>> future = std::async(std::launch::async, get_mean_chr_cov);
        futures.push_back(std::move(future));
    }

    // Loop through the futures and get the results
    for (auto& future : futures)
    {
        future.wait();
        std::tuple<uint32_t, uint32_t, std::unordered_map<uint32_t, int>> result = std::move(future.get());

        // Update the position count, cumulative depth, and merge the position-depth maps
        pos_count += std::get<0>(result);
        cum_depth += std::get<1>(result);
        this->mergePosDepthMaps(this->pos_depth_map, std::get<2>(result));
    }
    double mean_chr_cov = (double) cum_depth / (double) pos_count;

    return mean_chr_cov;
}

void CNVCaller::mergePosDepthMaps(std::unordered_map<uint32_t, int>& main_map, std::unordered_map<uint32_t, int>& map_update)
{
    // Merge the second depth map into the first
    main_map.reserve(main_map.size() + map_update.size());
    for (auto& pos_depth : map_update)
    {
        main_map[pos_depth.first] = std::move(pos_depth.second);
    }
}

double CNVCaller::calculateLog2Ratio(uint32_t start_pos, uint32_t end_pos, std::unordered_map<uint32_t, int> &pos_depth_map, double mean_chr_cov)
{
    // Use the position and depth map to calculate the log2 ratio
    double cum_depth = 0;
    int pos_count = 0;
    for (uint32_t i = start_pos; i <= end_pos; i++)
    {
        // Check if the position is in the map
        auto it = pos_depth_map.find(i);
        if (it == pos_depth_map.end())
        {
            continue;
        }
        int depth = pos_depth_map[i];
        pos_count++;
        cum_depth += depth;
    }

    // Calculate the window coverage log2 ratio (0 if no positions)
    double window_mean_cov = 0;
    if (pos_count > 0)
    {
        window_mean_cov = (double) cum_depth / (double) pos_count;
    }

    // Calculate the log2 ratio for the window
    // Avoid log2(0) by using a small value
    if (window_mean_cov == 0)
    {
        window_mean_cov = 0.0001;
    }
    double window_log2_ratio = log2(window_mean_cov / mean_chr_cov);

    return window_log2_ratio;
}

void CNVCaller::readSNPAlleleFrequencies(std::string chr, std::string filepath, SNPInfo& snp_info)
{
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
    pclose(index_fp);  // Close the process

    // Filter variants by depth, quality, and region
    if (this->input_data->getVerbose()) {
        std::cout << "Filtering SNPs by depth, quality, and region..." << std::endl;
    }

    // Check if a region was specified by the user
    std::string region_str = chr;
    if (this->input_data->isRegionSet())
    {
        std::pair<int32_t, int32_t> region = this->input_data->getRegion();
        region_str = chr + ":" + std::to_string(region.first) + "-" + std::to_string(region.second);
    }

    std::string filtered_snp_vcf_filepath = this->input_data->getOutputDir() + "/filtered_snps.vcf";
    std::string cmd = "bcftools view -r " + region_str + " -v snps -i 'QUAL > 30 && DP > 10 && FILTER = \"PASS\"' " + filepath + " > " + filtered_snp_vcf_filepath;
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

    pclose(fp);  // Close the process

    if (this->input_data->getVerbose()) {
        std::cout << "Finished extracting B-allele frequency data from filtered SNPs" << std::endl;
    }
}

void CNVCaller::getSNPPopulationFrequencies(std::string chr, SNPInfo& snp_info)
{
    // Get the population frequency file for the chromosome
    std::string pfb_filepath = this->input_data->getAlleleFreqFilepath(chr);
    if (pfb_filepath == "")
    {
        std::cout << "No population frequency file provided for chromosome " << chr << std::endl;
        return;
    }

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
    std::cout << "Reading population frequencies for chromosome " << chr << " from " << pfb_filepath << std::endl;

    // Get the start and end SNP positions for the chromosome (1-based
    // index)
    std::pair<int64_t, int64_t> snp_range = snp_info.getSNPRange(chr);
    int64_t snp_start = snp_range.first;
    int64_t snp_end = snp_range.second;
    if (this->input_data->isRegionSet())
    {
        // Get the user-defined region
        std::pair<int32_t, int32_t> region = this->input_data->getRegion();
        if (snp_start < region.first) {
            snp_start = region.first;
        } else if (snp_end > region.second) {
            snp_end = region.second;
        }
    }

    // Split region into chunks and get the population frequencies in parallel
    std::cout << "SNP range for chromosome " << chr << ": " << snp_start << "-" << snp_end << std::endl;
    int num_threads = this->input_data->getThreadCount();
    std::vector<std::string> region_chunks = splitRegionIntoChunks(chr_gnomad, snp_start, snp_end, num_threads);
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
            std::string filter_criteria = "INFO/variant_type=\"snv\" && " + AF_key + " >= " + std::to_string(MIN_PFB) + " && " + AF_key + " <= " + std::to_string(MAX_PFB);
            std::string cmd = \
                "bcftools query -r " + region_chunk + " -f '%POS\t%" + AF_key + "\n' -i '" + filter_criteria + "' " + pfb_filepath + " 2>/dev/null";

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
                    pos_pfb_map[pos] = pfb;  // Add the position and population frequency to the map
                }
            }
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
        future.wait();
        std::unordered_map<int, double> result = std::move(future.get());

        // Loop through the result and add to SNPInfo
        // printMessage("Adding population frequencies to SNPInfo...");
        for (auto& pair : result)
        {
            int pos = pair.first;
            double pfb = pair.second;

            // Add the population frequency to the SNPInfo
            this->snp_data_mtx.lock();
            snp_info.insertSNPPopulationFrequency(chr_snp, pos, pfb);
            this->snp_data_mtx.unlock();
            pfb_count++;

            // [TEST] Print 15 values
            if (pfb_count < 15)
            {
                printMessage("Population frequency for " + chr + ":" + std::to_string(pos) + " = " + std::to_string(pfb));
            }
        }
    }
}

void CNVCaller::saveSVCopyNumberToTSV(SNPData& snp_data, std::string filepath, std::string chr, int64_t start, int64_t end, std::string sv_type, double likelihood)
{
    // Open the TSV file for writing
    std::ofstream tsv_file(filepath);
    if (!tsv_file.is_open())
    {
        std::cerr << "ERROR: Could not open TSV file for writing: " << filepath << std::endl;
        exit(1);
    }

    // Ensure all values are valid, and print an error message if not
    if (chr == "" || start == 0 || end == 0 || sv_type == "")
    {
        std::cerr << "ERROR: Invalid SV information for TSV file: " << chr << ":" << start << "-" << end << " " << sv_type << std::endl;
        exit(1);
    }

    if (snp_data.pos.size() == 0)
    {
        std::cerr << "ERROR: No SNP data available for TSV file: " << chr << ":" << start << "-" << end << " " << sv_type << std::endl;
        exit(1);
    }

    if (likelihood == 0.0)
    {
        std::cerr << "ERROR: Invalid likelihood value for TSV file: " << chr << ":" << start << "-" << end << " " << sv_type << std::endl;
        exit(1);
    }

    if (sv_type.empty())
    {
        std::cerr << "ERROR: Invalid SV type for TSV file: " << chr << ":" << start << "-" << end << " " << sv_type << std::endl;
        exit(1);
    }

    // Format the size string in kb
    std::stringstream ss;
    ss << std::fixed << std::setprecision(6) << likelihood;
    std::string likelihood_str = ss.str();

    // Print SV information to the TSV file
    tsv_file << "SVTYPE=" << sv_type << std::endl;
    tsv_file << "POS=" << chr << ":" << start << "-" << end << std::endl;
    tsv_file << "HMM_LOGLH=" << likelihood_str << std::endl;

    // Write the header
    tsv_file << "chromosome\tposition\tsnp\tb_allele_freq\tlog2_ratio\tcnv_state\tpopulation_freq" << std::endl;

    // Write the data
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

        // If the SNP is not a SNP, then set the BAF to 0.0
        if (!is_snp)
        {
            baf = 0.0;
        }

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
