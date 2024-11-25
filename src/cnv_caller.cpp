
#include "cnv_caller.h"

#include <htslib/sam.h>

#include <htslib/vcf.h>
#include <htslib/hts.h>
#include <htslib/synced_bcf_reader.h>

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

CNVCaller::CNVCaller(InputData &input_data)
    : input_data(input_data)  // Initialize the input data
{
}

// Function to call the Viterbi algorithm for the CHMM
std::pair<std::vector<int>, double> CNVCaller::runViterbi(CHMM hmm, SNPData& snp_data)
{
    int data_count = (int) snp_data.pos.size();
    if (data_count == 0)
    {
        throw std::runtime_error("Error: No SNP data found for Viterbi algorithm.");
    }
    // std::lock_guard<std::mutex> lock(this->hmm_mtx);  // Lock the mutex for the HMM
    std::pair<std::vector<int>, double> state_sequence = testVit_CHMM(hmm, data_count, snp_data.log2_cov, snp_data.baf, snp_data.pfb);
    return state_sequence;
}

// Function to obtain SNP information for a region
std::pair<SNPData, bool> CNVCaller::querySNPRegion(std::string chr, uint32_t start_pos, uint32_t end_pos, SNPInfo& snp_info, std::unordered_map<uint32_t, int>& pos_depth_map, double mean_chr_cov)
{
    SNPData snp_data;
    bool snps_found = false;
    uint32_t window_size = (uint32_t)this->input_data.getWindowSize();

    // Query the SNPs for the entire region
    std::set<uint32_t> snp_pos;
    std::unordered_map<uint32_t, double> snp_baf;
    std::unordered_map<uint32_t, double> snp_pfb;
    this->querySNPs(chr, start_pos, end_pos, snp_pos, snp_baf, snp_pfb);
    // std::pair<std::vector<uint32_t>, std::vector<double>, std::vector<double>> snp_query = this->querySNPs(chr, start_pos, end_pos, snp_pos, snp_baf, snp_pfb);
    // std::vector<uint32_t>& snp_pos = std::get<0>(snp_query);
    // std::vector<double>& snp_pfb = std::get<1>(snp_query);
    // std::vector<double>& snp_baf = std::get<2>(snp_query);

    // Loop through the range of the SV region and query the SNPs in a sliding
    // window, then calculate the log2 ratio for each window
    for (uint32_t i = start_pos; i <= end_pos; i += window_size)
    {
        // std::cout << "Querying SNP region for " << chr << ":" << i << "-" << std::min(i + window_size - 1, end_pos) << std::endl;
        // Run a sliding non-overlapping window of size window_size across
        // the SV region and calculate the log2 ratio for each window
        uint32_t window_start = i;
        uint32_t window_end = std::min(i + window_size - 1, end_pos);

        // Get the SNP info for the window
        std::vector<uint32_t> snp_window_pos;
        std::vector<double> snp_window_bafs;
        std::vector<double> snp_window_pfbs;
        auto it_start = snp_pos.lower_bound(window_start);
        auto it_end = snp_pos.upper_bound(window_end);
        for (auto it = it_start; it != it_end; it++)
        {
            snp_window_pos.push_back(*it);
            snp_window_bafs.push_back(snp_baf[*it]);
            snp_window_pfbs.push_back(snp_pfb[*it]);
        }

        // Loop though the SNP positions and calculate the log2 ratio for
        // the window up to the SNP, then calculate the log2 ratio centered
        // at the SNP, and finally calculate the log2 ratio for the window
        // after the SNP, and continue until the end of the window
        // (If there are no SNPs in the window, then use the default BAF and
        // PFB values, and the coverage log2 ratio)

        // If no SNPs, then calculate the log2 ratio for the window
        if (snp_window_pos.size() == 0)
        {
            double window_log2_ratio = calculateLog2Ratio(window_start, window_end, pos_depth_map, mean_chr_cov);
            double pfb_default = 0.5;
            double baf_default = -1.0;  // Use -1.0 to indicate no BAF data
            this->updateSNPData(snp_data, (window_start + window_end) / 2, pfb_default, baf_default, window_log2_ratio, false);

        } else {
            snps_found = true;

            // Loop through the SNPs and calculate the log2 ratios
            // uint32_t bin_start = window_start;
            // uint32_t bin_end = 0;
            for (int j = 0; j < (int) snp_window_pos.size(); j++)
            {
                // Just use a window centered at the SNP position
                uint32_t bin_start = snp_window_pos[j] - window_size / 2;
                uint32_t bin_end = snp_window_pos[j] + window_size / 2;

                // Trim the bin start and end to 1/2 the distance from the
                // neighboring SNPs (or the start/end of the window)
                if (j > 0)
                {
                    bin_start = std::max(bin_start, (snp_window_pos[j-1] + snp_window_pos[j]) / 2);
                }

                if (j < (int) snp_window_pos.size() - 1)
                {
                    bin_end = std::min(bin_end, (snp_window_pos[j] + snp_window_pos[j+1]) / 2);
                }
                // std::cout << "bin_start: " << bin_start << std::endl;
                // std::cout << "bin_end: " << bin_end << std::endl;

                // Calculate the log2 ratio for the SNP bin
                double bin_cov = calculateLog2Ratio(bin_start, bin_end, pos_depth_map, mean_chr_cov);
                this->updateSNPData(snp_data, snp_window_pos[j], snp_window_pfbs[j], snp_window_bafs[j], bin_cov, true);
                // this->updateSNPData(snp_data, snp_pos, snp_window_pfbs[j], snp_window_bafs[j], bin_cov, true);

                // Update the previous bin start
                bin_start = bin_end + 1;
            }
        }
    }

    return std::make_pair(snp_data, snps_found);
}

std::tuple<double, SVType, std::string, bool> CNVCaller::runCopyNumberPrediction(std::string chr, const SVCandidate& candidate)
{
     // Get the start and end positions of the SV call
    uint32_t start_pos = std::get<0>(candidate);
    uint32_t end_pos = std::get<1>(candidate);

    // Run the Viterbi algorithm on SNPs in the SV region +/- 1/2
    // the SV length
    uint32_t sv_half_length = (end_pos - start_pos) / 2.0;
    // uint32_t snp_start_pos = std::max((uint32_t)1, start_pos - sv_length);
    // Prevent underflow (start_pos - sv_length) if start_pos < sv_length
    uint32_t snp_start_pos = start_pos > sv_half_length ? start_pos - sv_half_length : 1;
    uint32_t snp_end_pos = end_pos + sv_half_length;
    // std::cout << "CNP for " << chr << ":" << start_pos << "-" << end_pos << "(" << snp_start_pos << ", " << snp_end_pos << ")" << std::endl;
    // printMessage("Running copy number prediction for SV candidate " + chr + ":" + std::to_string(start_pos) + "-" + std::to_string(end_pos) + " with SNP region " + chr + ":" + std::to_string(snp_start_pos) + "-" + std::to_string(snp_end_pos) + "...");

    // Query the SNP region for the SV candidate
    std::pair<SNPData, bool> snp_call = querySNPRegion(chr, snp_start_pos, snp_end_pos, this->snp_info, this->pos_depth_map, this->mean_chr_cov);
    SNPData& sv_snps = snp_call.first;
    bool sv_snps_found = snp_call.second;

	/*
    if (sv_snps.pos.size() == 0) {
    	std::cerr << "ERROR [2]: No windows for SV " << chr << ":" << std::to_string((int)start_pos) << "-" << std::to_string((int)end_pos) << " (" << snp_start_pos << "," << snp_end_pos << std::endl;
    	continue;
    }
    */

    // Run the Viterbi algorithm
    // printMessage("[TEST] Running Viterbi algorithm for SV candidate " + chr + ":" + std::to_string(start_pos) + "-" + std::to_string(end_pos) + "...");
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

    // Combine counts for states 1 and 2, states 3 and 4, and states 5 and 6
    for (int i = 0; i < 6; i += 2)
    {
        // Combine counts for states 1 and 2, states 3 and 4, and states 5 and 6
        int state_count = std::count(sv_states.begin(), sv_states.end(), i+1) + std::count(sv_states.begin(), sv_states.end(), i+2);
        if (state_count > max_count)
        {
            max_state = i+1;  // Set the state to the first state in the pair (sequence remains intact)
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
    if (this->input_data.getSaveCNVData() && copy_number_change && (end_pos - start_pos) > 10000)
    {
        std::string cnv_type_str = getSVTypeString(predicted_cnv_type);
        std::string sv_filename = this->input_data.getOutputDir() + "/" + cnv_type_str + "_" + chr + "_" + std::to_string((int) start_pos) + "-" + std::to_string((int) end_pos) + "_SPLITALN.tsv";
        std::cout << "Saving SV split-alignment copy number predictions to " << sv_filename << std::endl;
        this->saveSVCopyNumberToTSV(sv_snps, sv_filename, chr, start_pos, end_pos, cnv_type_str, likelihood);
    }
    
    return std::make_tuple(likelihood, predicted_cnv_type, genotype, sv_snps_found);
}


void CNVCaller::runCIGARCopyNumberPrediction(std::string chr, std::set<SVCall> &sv_candidates, int min_length)
{
    CHMM& hmm = this->hmm;
    int window_size = this->input_data.getWindowSize();
    double mean_chr_cov = this->mean_chr_cov;  
    printMessage("Predicting CIGAR string copy number states for chromosome " + chr + "...");

    // Create a map with counts for each CNV type
    // std::map<int, int> cnv_type_counts;
    // for (int i = 0; i < 6; i++)
    // {
    //     cnv_type_counts[i] = 0;
    // }

    runCIGARCopyNumberPredictionChunk(chr, sv_candidates, hmm, window_size, mean_chr_cov);
    // // Split the SV candidates into chunks for each thread
    // int chunk_count = this->input_data.getThreadCount();
    // // std::vector<std::vector<SVCandidate>> sv_chunks = splitSVCandidatesIntoChunks(sv_candidates, chunk_count);
    // std::vector<std::set<SVCall>> sv_chunks = splitSVsIntoChunks(sv_candidates, chunk_count);

    // // Loop through each SV chunk and run the copy number prediction in parallel
    // // std::vector<std::future<SNPData>> futures;
    // std::vector<std::future<void>> futures;
    // for (auto& sv_chunk : sv_chunks)
    // {
    //     // Run the copy number prediction for the SV chunk
    //     futures.emplace_back(std::async(std::launch::async, &CNVCaller::runCIGARCopyNumberPredictionChunk, this, chr, std::ref(sv_chunk), std::ref(this->snp_info), hmm, window_size, mean_chr_cov, std::ref(this->pos_depth_map)));
    //     // futures.emplace_back(std::async(std::launch::async, &CNVCaller::runCIGARCopyNumberPredictionChunk, this, chr, std::ref(sv_chunk), std::ref(this->snp_info), hmm, window_size, mean_chr_cov, std::ref(this->pos_depth_map)));
    //     // std::async(std::launch::async, &CNVCaller::runCIGARCopyNumberPredictionChunk, this, chr, sv_chunk, std::ref(this->snp_info), hmm, window_size, mean_chr_cov, std::ref(this->pos_depth_map));
    // }

    // // Wait for all the futures to finish
    // int current_chunk = 0;
    // for (auto& future : futures)
    // {
    //     current_chunk++;
    //     try {
    //         future.wait();
    //         // SNPData chunk_snp_data = std::move(future.get());
    //         if (this->input_data.getVerbose())
    //         {
    //             printMessage("Finished processing SV chunk " + std::to_string(current_chunk) + " of " + std::to_string(chunk_count) + "...");
    //         }
    //     } catch (const std::exception& e) {
    //         printError("Error processing SV chunk " + std::to_string(current_chunk) + " of " + std::to_string(chunk_count) + ": " + e.what());
    //     }
    // }

    printMessage("Finished predicting copy number states for chromosome " + chr + "...");
}

void CNVCaller::runCIGARCopyNumberPredictionChunk(std::string chr, std::set<SVCall>& sv_chunk, CHMM hmm, int window_size, double mean_chr_cov)
{
    // printMessage("Running copy number prediction for " + std::to_string(sv_chunk.size()) + " SV candidates on chromosome " + chr + "...");
    // Map with counts for each CNV type
    std::map<int, int> cnv_type_counts;
    for (int i = 0; i < 6; i++)
    {
        cnv_type_counts[i] = 0;
    }
    
    // Loop through each SV candidate and predict the copy number state
    for (auto& sv_call : sv_chunk)
    {

        // Get the SV candidate
        // const SVCandidate& candidate = sv_call;
        // int64_t start_pos = std::get<0>(candidate);
        // int64_t end_pos = std::get<1>(candidate);
        uint32_t start_pos = sv_call.start;
        uint32_t end_pos = sv_call.end;
        
        // Error if start > end
        if (start_pos >= end_pos)
        {
        	std::cerr << "Position error for CIGAR SV at " << chr << ":" << start_pos << "-" << end_pos << std::endl;
        	continue;
        }

        // Skip if not the minimum length for CNV predictions
        if ((end_pos - start_pos) < (uint32_t)this->input_data.getMinCNVLength())
        {
            continue;
        }

    	// std::cout << "CIGAR SV at " << chr << ":" << start_pos << "-" << end_pos << std::endl;

        // Get the depth at the start position. This is used as the FORMAT/DP
        // value in the VCF file
        // int dp_value = pos_depth_map[start_pos];
        // this->updateDPValue(sv_candidates, sv_call, dp_value);

        // Loop through the SV region +/- 1/2 SV length and run copy number
        // predictions
        uint32_t sv_half_length = (end_pos - start_pos) / 2.0;
        uint32_t snp_start_pos = start_pos > sv_half_length ? start_pos - sv_half_length : 1;
        uint32_t snp_end_pos = end_pos + sv_half_length;
        // std::cout << "CIGAR sv_half_length:" << sv_half_length << std::endl;
        // std::cout << "CIGAR SV query at " << chr << ":" << query_start << "-" << query_end << std::endl;

        // printMessage("Querying SNPs for SV " + chr + ":" + std::to_string((int)start_pos) + "-" + std::to_string((int)end_pos) + ", qstart = " + std::to_string(query_start) + ", qend = " + std::to_string(query_end));
        std::pair<SNPData, bool> snp_call = this->querySNPRegion(chr, snp_start_pos, snp_end_pos, snp_info, this->pos_depth_map, mean_chr_cov);
        // printMessage("Finished querying SNPs for SV " + chr + ":" + std::to_string((int)start_pos) + "-" + std::to_string((int)end_pos));
        SNPData& sv_snps = snp_call.first;
        bool snps_found = snp_call.second;

        // Run the Viterbi algorithm
        // printMessage("[TEST2] Running Viterbi algorithm for SV " + chr + ":" + std::to_string((int)start_pos) + "-" + std::to_string((int)end_pos) + "...");
        
        if (sv_snps.pos.size() == 0) {
        	std::cerr << "ERROR: No windows for SV " << chr << ":" << start_pos << "-" << end_pos << " (" << snp_start_pos << "," << snp_end_pos << std::endl;
        	continue;
        }
        
        std::pair<std::vector<int>, double> prediction = runViterbi(hmm, sv_snps);
        std::vector<int>& state_sequence = prediction.first;
        double likelihood = prediction.second;
        // printMessage("Finished running Viterbi algorithm for SV " + chr + ":" + std::to_string((int)start_pos) + "-" + std::to_string((int)end_pos) + "...");

        // Get all the states in the SV region
        // printMessage("Getting states for SV " + chr + ":" + std::to_string((int)start_pos) + "-" + std::to_string((int)end_pos) + "...");
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
        SVType updated_sv_type = getSVTypeFromCNState(max_state);
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

        // Update the SV copy number data if not unknown
        // printMessage("Updating SV copy number data for SV " + chr + ":" + std::to_string((int)start_pos) + "-" + std::to_string((int)end_pos) + "...");
        if (updated_sv_type != SVType::UNKNOWN)
        {
            std::string sv_type_str = getSVTypeString(updated_sv_type);
            addSVCall(sv_chunk, sv_call.start, sv_call.end, sv_type_str, ".", data_type, genotype, likelihood);
            // std::string sv_type_str = getSVTypeString(updated_sv_type);
            // sv_call.sv_type = sv_type_str;
            // sv_call.data_type += "," + data_type;
            // sv_call.genotype = genotype;
            // sv_call.hmm_likelihood = likelihood;
        }
        // this->updateSVCopyNumber(sv_candidates, sv_call, cnv_type, data_type, genotype, likelihood);

        // Save the SV calls as a TSV file if enabled, if the SV type is
        // known, and the length is greater than 10 kb
        // SVType updated_sv_type = sv_candidates[sv_call].sv_type;
        if (this->input_data.getSaveCNVData() && updated_sv_type != SVType::UNKNOWN && (end_pos - start_pos) > 10000)
        {
            // Add the state sequence to the SNP data (avoid copying the data)
            sv_snps.state_sequence = std::move(state_sequence);

            // Save the SV calls as a TSV file
            std::string cnv_type_str = getSVTypeString(updated_sv_type);
            std::string sv_filename = this->input_data.getOutputDir() + "/" + cnv_type_str + "_" + chr + "_" + std::to_string((int) start_pos) + "-" + std::to_string((int) end_pos) + "_CIGAR.tsv";
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

std::vector<std::string> CNVCaller::splitRegionIntoChunks(std::string chr, uint32_t start_pos, uint32_t end_pos, int chunk_count)
{
    // Split the region into chunks
    std::vector<std::string> region_chunks;
    uint32_t region_length = end_pos - start_pos + 1;
    uint32_t chunk_size = std::ceil((double) region_length / (double) chunk_count);
    uint32_t chunk_start = start_pos;
    uint32_t chunk_end = 0;
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

void CNVCaller::loadChromosomeData(std::string chr)
{
    std::string hmm_filepath = this->input_data.getHMMFilepath();
    std::cout << "Reading HMM from file: " << hmm_filepath << std::endl;
    this->hmm = ReadCHMM(hmm_filepath.c_str());

    printMessage("Calculating mean chromosome coverage for " + chr + "...");
    this->mean_chr_cov = calculateMeanChromosomeCoverage(chr);
    //this->mean_chr_cov = 30.0;
    printMessage("Mean chromosome coverage for " + chr + ": " + std::to_string(mean_chr_cov));

    // std::cout << "Reading SNP allele frequencies for chromosome " << chr << " from VCF file..." << std::endl;
    // std::string snp_filepath = this->input_data.getSNPFilepath();
    // readSNPAlleleFrequencies(chr, snp_filepath, this->snp_info);

    // std::cout << "Obtaining SNP population frequencies for chromosome " << chr << "..." << std::endl;
    // getSNPPopulationFrequencies(chr, this->snp_info);
    // std::cout << "Finished loading chromosome data for " << chr << std::endl;
}

// Calculate the mean chromosome coverage
double CNVCaller::calculateMeanChromosomeCoverage(std::string chr)
{
    // Open the BAM file
    std::string bam_filepath = this->input_data.getShortReadBam();
    samFile *bam_file = sam_open(bam_filepath.c_str(), "r");
    if (!bam_file)
    {
        throw std::runtime_error("ERROR: Could not open BAM file: " + bam_filepath);
    }

    // Enable multi-threading
    hts_set_threads(bam_file, this->input_data.getThreadCount());

    // Read the header
    bam_hdr_t *bam_header = sam_hdr_read(bam_file);
    if (!bam_header)
    {
        sam_close(bam_file);
        throw std::runtime_error("ERROR: Could not read header from BAM file: " + bam_filepath);
    }

    // Load the index
    hts_idx_t *bam_index = sam_index_load(bam_file, bam_filepath.c_str());
    if (!bam_index)
    {
        bam_hdr_destroy(bam_header);
        sam_close(bam_file);
        throw std::runtime_error("ERROR: Could not load index for BAM file: " + bam_filepath);
    }

    // Create an iterator for the chromosome
    hts_itr_t *bam_iter = sam_itr_querys(bam_index, bam_header, chr.c_str());
    if (!bam_iter)
    {
        hts_idx_destroy(bam_index);
        bam_hdr_destroy(bam_header);
        sam_close(bam_file);
        throw std::runtime_error("ERROR: Could not create iterator for chromosome: " + chr);
    }

    // Initialize the record
    bam1_t *bam_record = bam_init1();
    if (!bam_record)
    {
        hts_itr_destroy(bam_iter);
        hts_idx_destroy(bam_index);
        bam_hdr_destroy(bam_header);
        sam_close(bam_file);
        throw std::runtime_error("ERROR: Could not initialize BAM record.");
    }

    // Iterate through the chromosome and update the depth map
    std::unordered_map<uint32_t, int> chr_pos_depth_map;
    while (sam_itr_next(bam_file, bam_iter, bam_record) >= 0)
    {
        // Ignore UNMAP, SECONDARY, QCFAIL, and DUP reads
        if (bam_record->core.flag & BAM_FUNMAP || bam_record->core.flag & BAM_FSECONDARY || bam_record->core.flag & BAM_FQCFAIL || bam_record->core.flag & BAM_FDUP)
        {
            continue;
        }
        
        // Parse the CIGAR string to get the depth (match, sequence match, and
        // mismatch)
        // uint32_t depth = 0;
        uint32_t pos = bam_record->core.pos + 1;  // 0-based to 1-based
        uint32_t ref_pos = pos;
        uint32_t cigar_len = bam_record->core.n_cigar;
        uint32_t *cigar = bam_get_cigar(bam_record);
        for (uint32_t i = 0; i < cigar_len; i++)
        {
            uint32_t op = bam_cigar_op(cigar[i]);
            uint32_t op_len = bam_cigar_oplen(cigar[i]);
            if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF)
            {
                // Update the depth for each position in the alignment
                for (uint32_t j = 0; j < op_len; j++)
                {
                    chr_pos_depth_map[ref_pos + j]++;
                }
            }
            
            // Update the reference coordinate based on the CIGAR operation
            // https://samtools.github.io/hts-specs/SAMv1.pdf
            if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP || op == BAM_CEQUAL || op == BAM_CDIFF) {
                ref_pos += op_len;
            } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP || op == BAM_CPAD) {
                // Do nothing
            } else {
                throw std::runtime_error("ERROR: Unknown CIGAR operation: " + std::to_string(op));
            }
        }
    }

    // Clean up
    bam_destroy1(bam_record);
    hts_itr_destroy(bam_iter);
    hts_idx_destroy(bam_index);
    bam_hdr_destroy(bam_header);
    sam_close(bam_file);

    // Calculate the mean chromosome coverage for positions with non-zero depth
    uint64_t cum_depth = 0;
    uint32_t pos_count = 0;
    for (auto& pos_depth : chr_pos_depth_map)
    {
        // if (pos_depth.second > 0)
        // {
        //     cum_depth += pos_depth.second;
        //     pos_count++;
        // } else {
        //     std::cout << "Zero depth at position " << pos_depth.first << std::endl;
        // }
        cum_depth += pos_depth.second;
        pos_count++;
    }

    double mean_chr_cov = (double) cum_depth / (double) pos_count;

    // Update the position depth map
    this->pos_depth_map = std::move(chr_pos_depth_map);

    return mean_chr_cov;
}

void CNVCaller::mergePosDepthMaps(std::unordered_map<uint32_t, int>& main_map, std::unordered_map<uint32_t, int>& map_update)
{
    // Merge the second depth map into the first
    for (auto& pos_depth : map_update)
    {
        main_map[pos_depth.first] = pos_depth.second;
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

void CNVCaller::readSNPAlleleFrequencies(std::string chr, uint32_t start_pos, uint32_t end_pos, std::set<uint32_t>& snp_pos, std::unordered_map<uint32_t, double>& snp_baf)
{
    // Get the SNP file path
    std::string snp_filepath = this->input_data.getSNPFilepath();
    if (snp_filepath.empty())
    {
        throw std::runtime_error("ERROR: SNP file path is empty.");
    }

    // Initialize the synced reader
    bcf_srs_t *snp_reader = bcf_sr_init();
    if (!snp_reader)
    {
        throw std::runtime_error("ERROR: Could not initialize SNP reader.");
    }

    // Set the region
    std::string region_str = chr + ":" + std::to_string(start_pos) + "-" + std::to_string(end_pos);
    if (bcf_sr_set_regions(snp_reader, region_str.c_str(), 0) < 0)
    {
        bcf_sr_destroy(snp_reader);
        throw std::runtime_error("ERROR: Could not set region for SNP reader: " + region_str);
    }

    // Set multi-threading
    int thread_count = this->input_data.getThreadCount();
    bcf_sr_set_threads(snp_reader, thread_count);

    // Enable index usage
    snp_reader->require_index = 1;

    // Add the SNP file to the reader
    if (bcf_sr_add_reader(snp_reader, snp_filepath.c_str()) < 0)
    {
        bcf_sr_destroy(snp_reader);
        throw std::runtime_error("ERROR: Could not add SNP file to reader: " + snp_filepath);
    }

    // Get the header
    bcf_hdr_t *snp_header = bcf_sr_get_header(snp_reader, 0);
    if (!snp_header)
    {
        bcf_sr_destroy(snp_reader);
        throw std::runtime_error("ERROR: Could not get header for SNP reader.");
    }

    // std::cout << "Iterating through SNPs in region " << region_str << "..." << std::endl;
    int print_count = 0;
    int record_count = 0;
    int duplicate_count = 0;
    uint32_t last_pos = 0;
    while (bcf_sr_next_line(snp_reader) > 0)
    {
        if (!bcf_sr_has_line(snp_reader, 0))
        {
            continue;
        }
        bcf1_t *snp_record = bcf_sr_get_line(snp_reader, 0);
        if (snp_record)
        {
            record_count++;
            uint32_t pos = (uint32_t)snp_record->pos + 1;

            // Skip if 3 or more duplicate positions found
            // if (pos == last_pos)
            // {
            //     duplicate_count++;
            //     if (duplicate_count >= 10)
            //     {
            //         std::cerr << "ERROR: 3 or more duplicate positions found in SNP file at " << chr << ":" << pos << std::endl;
            //         break;
            //     }
            // } else {
            //     duplicate_count = 0;
            // }
            // last_pos = pos;

            // Skip if not a SNP
            if (!bcf_is_snp(snp_record))
            {
                continue;
            }

            // Get the QUAL, DP, and AD values
            float qual = snp_record->qual;
            if (bcf_float_is_missing(qual))
            {
                // std::cerr << "ERROR: QUAL value is missing for SNP at " << chr << ":" << pos << std::endl;
            }
            // Skip if quality is less than 30
            if (qual <= 30)
            {
                continue;
            }

            // Extract DP from FORMAT field
            int32_t *dp = 0;
            int dp_count = 0;
            int dp_ret = bcf_get_format_int32(snp_header, snp_record, "DP", &dp, &dp_count);
            bool dp_skip = false;
            if (dp_ret < 0)
            {
                // std::cerr << "ERROR: Could not get DP value for SNP at " << chr << ":" << pos << std::endl;
            } else {
                // Skip if depth is not greater than 10
                for (int i = 0; i < dp_count; i++)
                {
                    if (dp[i] <= 10)
                    {
                        dp_skip = true;
                        break;
                    }
                }
            }
            free(dp);
            if (dp_skip)
            {
                continue;
            }

            // Skip if the SNP does not pass the filter
            if (bcf_has_filter(snp_header, snp_record, const_cast<char*>("PASS")) != 1)
            {
                continue;
            }

            // Extract AD from FORMAT field
            int32_t *ad = 0;
            int ad_count = 0;
            int ad_ret = bcf_get_format_int32(snp_header, snp_record, "AD", &ad, &ad_count);

            // Skip if AD value is missing
            if (ad_ret < 0)
            {
                // std::cerr << "ERROR: AD value is missing for SNP at " << chr
                // << ":" << pos << std::endl;
                throw std::runtime_error("ERROR: AD value is missing for SNP at " + chr + ":" + std::to_string(pos));
            }

            // Calculate the B-allele frequency (BAF)
            double baf = 0.0;
            double ad0 = 0.0;
            double ad1 = 0.0;
            for (int i = 0; i < ad_count; i++)
            {
                if (i == 0)
                {
                    ad0 = (double) ad[i];
                } else if (i == 1) {
                    ad1 = (double) ad[i];
                }
            }
            free(ad);
            baf = ad1 / (ad0 + ad1);

            // Insert the SNP position and BAF into the maps
            snp_pos.insert(pos);
            snp_baf[pos] = baf;

            // Print the SNP position and BAF
            // std::cout << "SNP: " << chr << ":" << pos << ", BAF: " << baf << "(REF=" << ad0 << ",ALT=" << ad1 << ")" << std::endl;
            // print_count++;
            // if (print_count < 10)
            // {
            //     std::cout << "SNP: " << chr << ":" << pos << ", BAF: " << baf << "(REF=" << ad0 << ",ALT=" << ad1 << ")" << std::endl;
            //     print_count++;
            // }
        }
    }

    // Clean up
    bcf_sr_destroy(snp_reader);

    // std::cout << "Opening SNP file: " << snp_filepath << std::endl;
    // htsFile *snp_file = bcf_open(snp_filepath.c_str(), "r");
    // if (!snp_file)
    // {
    //     throw std::runtime_error("ERROR: Could not open SNP file: " + snp_filepath);
    // }

    // // Enable multi-threading
    // hts_set_threads(snp_file, thread_count);

    // // Read the header
    // bcf_hdr_t *snp_header = bcf_hdr_read(snp_file);
    // if (!snp_header)
    // {
    //     bcf_close(snp_file);
    //     throw std::runtime_error("ERROR: Could not read header from SNP file: " + snp_filepath);
    // }

    // // Load the index
    // hts_idx_t *snp_index = bcf_index_load(snp_filepath.c_str());
    // if (!snp_index)
    // {
    //     bcf_hdr_destroy(snp_header);
    //     bcf_close(snp_file);
    //     throw std::runtime_error("ERROR: Could not load index for SNP file: " + snp_filepath);
    // }

    // // Construct the region string
    // std::string region_str = chr + ":" + std::to_string(start_pos) + "-" + std::to_string(end_pos);
    // hts_itr_t *snp_iter = bcf_itr_querys(snp_index, snp_header, region_str.c_str());
    // if (!snp_iter)
    // {
    //     hts_idx_destroy(snp_index);
    //     bcf_hdr_destroy(snp_header);
    //     bcf_close(snp_file);
    //     throw std::runtime_error("ERROR: Could not create iterator for SNP region: " + region_str);
    // }

    // // Set up the record
    // bcf1_t *snp_record = bcf_init();
    // if (!snp_record)
    // {
    //     bcf_hdr_destroy(snp_header);
    //     bcf_close(snp_file);
    //     throw std::runtime_error("ERROR: Could not initialize SNP record.");
    // }

    // // Read the SNPs in the chromosome region
    // int print_count = 0;
    // while (bcf_itr_next(snp_file, snp_iter, snp_record) >= 0)
    // {
    //     // Get the position and B-allele frequency (BAF) from the SNP record
    //     uint32_t pos = snp_record->pos + 1;  // 0-based to 1-based

    //     // Get QUAL, DP, and AD values
    //     float qual = snp_record->qual;
    //     if (bcf_float_is_missing(qual))
    //     {
    //         std::cerr << "ERROR: QUAL value is missing for SNP at " << chr << ":" << pos << std::endl;
    //     }
    //     // Skip if quality is less than 30
    //     if (qual <= 30)
    //     {
    //         continue;
    //     }

    //     // Get FILTER status
    //     int pass_id = bcf_hdr_id2int(snp_header, BCF_DT_ID, "PASS");
    //     if (pass_id == -1)
    //     {
    //         std::cerr << "ERROR: Could not get PASS ID for SNP at " << chr << ":" << pos << std::endl;
    //     }
    //     std::string pass_filter = "PASS";
    //     if (bcf_has_filter(snp_header, snp_record, const_cast<char*>(pass_filter.c_str())) != 1)
    //     {
    //         // Skip if the SNP does not pass the filter
    //         continue;
    //     }

    //     // Extract DP from INFO field
    //     int32_t dp = 0;
    //     int dp_count = bcf_get_info_int32(snp_header, snp_record, "DP", &dp, &dp_count);
    //     if (dp_count != 1)
    //     {
    //         std::cerr << "ERROR: Could not get DP value for SNP at " << chr << ":" << pos << std::endl;
    //     }
    //     // Skip if depth is not greater than 10
    //     if (dp <= 10)
    //     {
    //         continue;
    //     }

    //     // Skip if not a SNP
    //     if (!bcf_is_snp(snp_record))
    //     {
    //         continue;
    //     }

    //     // Extract AD from FORMAT field
    //     int32_t ad[2] = {0, 0};
    //     int ad_count = 0;
    //     int ad_ret = bcf_get_format_int32(snp_header, snp_record, "AD", &ad, &ad_count);
    //     // if (ad_count != 2)
    //     // {
    //     //     std::cerr << "ERROR: Could not get AD value for SNP at " << chr << ":" << pos << std::endl;
    //     // }

    //     // Calculate the BAF
    //     if (ad_ret > 0 && ad_count > 0)
    //     {
    //         double baf = (double) ad[1] / (double) (ad[0] + ad[1]);
    //         snp_pos.insert(pos);
    //         snp_baf[pos] = baf;

    //         // Print the SNP position and BAF
    //         if (print_count < 10)
    //         {
    //             std::cout << "SNP: " << chr << ":" << pos << ", BAF: " << baf << std::endl;
    //             print_count++;
    //         }
    //     }
    // }

    // // Clean up
    // bcf_destroy(snp_record);
    // hts_itr_destroy(snp_iter);
    // hts_idx_destroy(snp_index);
    // bcf_hdr_destroy(snp_header);
    // bcf_close(snp_file);

    // // Check that the SNP file is sorted by running bcftools index and reading
    // // the error output
    // std::string index_cmd = "bcftools index " + filepath + " 2>&1 | grep -i error";
    // if (this->input_data.getVerbose()) {
    //     std::cout << "Command: " << index_cmd << std::endl;
    // }

    // // Open a pipe to read the output of the command
    // FILE *index_fp = popen(index_cmd.c_str(), "r");
    // if (index_fp == NULL)
    // {
    //     std::cerr << "ERROR: Could not open pipe for command: " << index_cmd << std::endl;
    //     exit(1);
    // }

    // // Read the output of the command
    // const int error_size = 256;
    // char index_error[error_size];
    // while (fgets(index_error, error_size, index_fp) != NULL)
    // {
    //     std::cerr << "ERROR: " << index_error << std::endl;
    //     exit(1);
    // }
    // pclose(index_fp);  // Close the process

    // // Filter variants by depth, quality, and region
    // if (this->input_data.getVerbose()) {
    //     std::cout << "Filtering SNPs by depth, quality, and region..." << std::endl;
    // }

    // // Check if a region was specified by the user
    // std::string region_str = chr;
    // if (this->input_data.isRegionSet())
    // {
    //     std::pair<int32_t, int32_t> region = this->input_data.getRegion();
    //     region_str = chr + ":" + std::to_string(region.first) + "-" + std::to_string(region.second);
    // }

    // std::string filtered_snp_vcf_filepath = this->input_data.getOutputDir() + "/filtered_snps.vcf";
    // int thread_count = this->input_data.getThreadCount();
    // // std::string cmd = "bcftools view -r " + region_str + " -v snps -i 'QUAL > 30 && DP > 10 && FILTER = \"PASS\"' " + filepath + " > " + filtered_snp_vcf_filepath;
    // std::string cmd = "bcftools view --threads " + std::to_string(thread_count) + " -r " + region_str + " -v snps -i 'QUAL > 30 && DP > 10 && FILTER = \"PASS\"' " + filepath + " > " + filtered_snp_vcf_filepath;
    // if (this->input_data.getVerbose()) {
    //     std::cout << "Filtering SNPs by depth and quality..." << std::endl;
    //     std::cout << "Command: " << cmd << std::endl;
    // }
    // system(cmd.c_str());
    
    // if (this->input_data.getVerbose()) {
    //     std::cout << "Filtered SNPs written to " << filtered_snp_vcf_filepath << std::endl;
    // }

    // // Extract B-allele frequency data from the VCF file and sort by chromosome
    // // and position
    // if (this->input_data.getVerbose()) {
    //     std::cout << "Extracting B-allele frequency data from filtered SNPs..." << std::endl;
    // }
    // cmd = "bcftools query -f '%POS,[%AD]\n' " + filtered_snp_vcf_filepath + " 2>/dev/null";
    // FILE *fp = popen(cmd.c_str(), "r");
    // if (fp == NULL)
    // {
    //     std::cerr << "ERROR: Could not open pipe for command: " << cmd << std::endl;
    //     exit(1);
    // }

    // // Read the reference and alternate allele depths from the VCF file
    // std::string alt_allele = "";  // Alternate allele
    // uint32_t pos = 0;
    // int ref_ad = 0;
    // int alt_ad = 0;
    // const int line_size = 1024;
    // char line[line_size];  // Line buffer
    // std::vector<int64_t> locations;
    // std::vector<double> bafs;
    // std::string chr_no_prefix = removeChrPrefix(chr);
    // while (fgets(line, line_size, fp) != NULL)
    // {
    //     // Parse the line
    //     char *tok = strtok(line, ",");  // Tokenize the line
    //     int col = 0;  // Column index
    //     while (tok != NULL)
    //     {
    //         // Get the position from column 2
    //         if (col == 0)
    //         {
    //             pos = (uint32_t)atoi(tok);
    //         }

    //         // Get the AD for the reference allele from column 3
    //         else if (col == 1)
    //         {
    //             ref_ad = atoi(tok);
    //         }

    //         // Get the AD for the non-reference allele from column 4
    //         else if (col == 2)
    //         {
    //             alt_ad = atoi(tok);
    //         }

    //         // Move to the next token
    //         tok = strtok(NULL, ",");
    //         col++;
    //     }

    //     // Calculate the B-allele frequency (BAF) as the ratio of the alternate
    //     // allele depth to the total depth (reference + alternate)
    //     double baf = (double) alt_ad / (double) (ref_ad + alt_ad);

    //     // Add a new location and BAF value to the chromosome's SNP data
    //     // (population frequency and log2 ratio will be added later)
    //     // snp_info.insertSNPAlleleFrequency(chr_no_prefix, pos, baf);
    //     this->snp_baf_map[pos] = baf;
    //     this->snp_baf_keys.insert(pos);
    // }

    // pclose(fp);  // Close the process

    if (this->input_data.getVerbose()) {
        std::cout << "Finished extracting B-allele frequency data from filtered SNPs" << std::endl;
    }
}

// void CNVCaller::getSNPPopulationFrequencies(std::string chr, SNPInfo& snp_info)
// {
//     // Get the population frequency file for the chromosome
//     std::string pfb_filepath = this->input_data.getAlleleFreqFilepath(chr);
//     if (pfb_filepath.empty())
//     {
//         std::cout << "No population frequency file provided for chromosome " << chr << std::endl;
//         return;
//     }
    
//     // Determine the ethnicity-specific allele frequency key
//     std::string AF_key = "AF";
//     if (this->input_data.getEthnicity() != "")
//     {
//         AF_key += "_" + this->input_data.getEthnicity();
//     }

//     // Check if the filepath uses the 'chr' prefix notations based on the
//     // chromosome name (*.chr1.vcf.gz vs *.1.vcf.gz)
//     std::string chr_gnomad = chr;  // gnomAD data may or may not have the 'chr' prefix
//     std::string chr_prefix = "chr";
//     if (pfb_filepath.find(chr_prefix) == std::string::npos)
//     {
//         // Remove the 'chr' prefix from the chromosome name
//         if (chr_gnomad.find(chr_prefix) != std::string::npos)
//         {
//             chr_gnomad = chr_gnomad.substr(chr_prefix.length());
//         }
//     } else {
//         // Add the 'chr' prefix to the chromosome name
//         if (chr_gnomad.find(chr_prefix) == std::string::npos)
//         {
//             chr_gnomad = chr_prefix + chr;
//         }
//     }

//     // Remove the 'chr' prefix from the chromosome name for SNP data. All
//     // SNP data in this program does not use the 'chr' prefix
//     std::string chr_no_prefix = removeChrPrefix(chr);

//     std::cout << "Reading population frequencies for chromosome " << chr << " from " << pfb_filepath << std::endl;
//     int thread_count = this->input_data.getThreadCount();

//     // Open the population frequency file
//     std::cout << "Opening population frequency file: " << pfb_filepath << std::endl;
//     htsFile *pfb_file = hts_open(pfb_filepath.c_str(), "r");
//     if (!pfb_file)
//     {
//         throw std::runtime_error("ERROR: Could not open population frequency file: " + pfb_filepath);
//     }

//     // Enable multi-threading
//     std::cout << "Setting number of threads to " << thread_count << std::endl;
//     hts_set_threads(pfb_file, thread_count);

//     // Read the header
//     std::cout << "Reading header from population frequency file..." << std::endl;
//     bcf_hdr_t *pfb_header = bcf_hdr_read(pfb_file);
//     if (!pfb_header)
//     {
//         bcf_close(pfb_file);
//         throw std::runtime_error("ERROR: Could not read header from population frequency file: " + pfb_filepath);
//     }

//     // Set up the record
//     std::cout << "Initializing BCF record..." << std::endl;
//     bcf1_t *pfb_record = bcf_init();
//     if (!pfb_record)
//     {
//         bcf_hdr_destroy(pfb_header);
//         bcf_close(pfb_file);
//         throw std::runtime_error("ERROR: Could not initialize BCF record.");
//     }

//     // Read the population frequencies for the chromosome
//     std::cout << "[TEST] Reading population frequencies for chromosome " << chr << " (AF_key = " << AF_key << ")..." << std::endl;
//     int print_count = 0;
//     while (bcf_read(pfb_file, pfb_header, pfb_record) == 0)
//     {
//         // Get the chromosome and position
//         // std::cout << "Reading record..." << std::endl;
//         uint32_t pos = pfb_record->pos + 1;  // 0-based to 1-based

//         // Skip if not a SNP, or if the position is not in the BAF map
//         if (!bcf_is_snp(pfb_record) || this->snp_baf_keys.count(pos) == 0)
//         {
//             continue;
//         }

//         // Get the population frequency for the SNP
//         float *pfb_f = NULL;
//         int count = 0;
//         int pfb_status = bcf_get_info_float(pfb_header, pfb_record, AF_key.c_str(), &pfb_f, &count);
//         if (pfb_status < 0 || count == 0)
//         {
//             std::cout << "Field " << AF_key << " not found, or count is 0" << std::endl;
//             continue;
//         }
//         double pfb = (double) pfb_f[0];
//         free(pfb_f);

//         // Continue if the population frequency is outside the threshold
//         if (pfb <= MIN_PFB || pfb >= MAX_PFB)
//         {
//             continue;
//         }

//         // Add the population frequency to the SNP data
//         // snp_info.insertSNPPopulationFrequency(chr_no_prefix, pos, pfb);
//         if (this->snp_pfb_map.find(pos) == this->snp_pfb_map.end())
//         {
//             this->snp_pfb_map[pos] = pfb;
//         } else {
//             // Keep the larger population frequency
//             if (pfb > this->snp_pfb_map[pos])
//             {
//                 this->snp_pfb_map[pos] = pfb;
//             }
//         }
//         if (print_count < 10)
//         {
//             std::cout << "Population frequency for " << chr << ":" << pos << " = " << pfb << std::endl;
//             print_count++;
//         }
//     }
//     std::cout << "Finished reading population frequencies for chromosome " << chr << std::endl;
// }

void CNVCaller::readSNPPopulationFrequencies(std::string chr, uint32_t start_pos, uint32_t end_pos, std::unordered_map<uint32_t, double>& snp_pfb_map)
{
    // Get the population frequency file for the chromosome
    std::string pfb_filepath = this->input_data.getAlleleFreqFilepath(chr);
    if (pfb_filepath == "")
    {
        std::cout << "No population frequency file provided for chromosome " << chr << std::endl;
        return;
    }
    
    // Determine the ethnicity-specific allele frequency key
    std::string AF_key = "AF";
    if (this->input_data.getEthnicity() != "")
    {
        AF_key += "_" + this->input_data.getEthnicity();
    }

    // Check if the filepath uses the 'chr' prefix notations based on the
    // chromosome name (*.chr1.vcf.gz vs *.1.vcf.gz)
    std::string chr_gnomad = chr;  // gnomAD data may or may not have the 'chr' prefix
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

    // Remove the 'chr' prefix from the chromosome name for SNP data. All
    // SNP data in this program does not use the 'chr' prefix
    std::string chr_no_prefix = removeChrPrefix(chr);
    int thread_count = this->input_data.getThreadCount();

    // Initialize the synced reader
    bcf_srs_t *pfb_reader = bcf_sr_init();
    if (!pfb_reader)
    {
        throw std::runtime_error("ERROR: Could not initialize synced reader for population frequency file: " + pfb_filepath);
    }

    // Set the region for the synced reader
    std::string region_str = chr_gnomad + ":" + std::to_string(start_pos) + "-" + std::to_string(end_pos);
    if (bcf_sr_set_regions(pfb_reader, region_str.c_str(), 0) < 0)
    {
        bcf_sr_destroy(pfb_reader);
        throw std::runtime_error("ERROR: Could not set region for synced reader: " + region_str);
    }

    // Set multi-threading
    bcf_sr_set_threads(pfb_reader, thread_count);

    // Enable index usage
    pfb_reader->require_index = 1;

    // Add the population frequency file to the synced reader
    if (bcf_sr_add_reader(pfb_reader, pfb_filepath.c_str()) < 0)
    {
        bcf_sr_destroy(pfb_reader);
        throw std::runtime_error("ERROR: Could not add population frequency file to synced reader: " + pfb_filepath);
    }

    // Get the header
    bcf_hdr_t *pfb_header = bcf_sr_get_header(pfb_reader, 0);
    if (!pfb_header)
    {
        bcf_sr_destroy(pfb_reader);
        throw std::runtime_error("ERROR: Could not get header for population frequency file: " + pfb_filepath);
    }

    int test_count = 0;
    int record_count = 0;
    while (bcf_sr_next_line(pfb_reader) > 0)
    {
        if (!bcf_sr_has_line(pfb_reader, 0))
        {
            continue;
        }
        // std::cout << "Reading record..." << std::endl;
        // pfb_record = bcf_sr_get_line(pfb_reader, 0);
        bcf1_t *pfb_record = bcf_sr_get_line(pfb_reader, 0);
        // Do something with the record
        if (pfb_record)
        {
            record_count++;
            // Skip if not a SNP
            if (!bcf_is_snp(pfb_record))
            {
                // std::cout << "Skipping non-SNP at " << chr << ":" << pfb_record->pos << std::endl;
                continue;
            }

            uint32_t pos = (uint32_t) pfb_record->pos + 1;  // 0-based to 1-based

            // Get the population frequency for the SNP
            float *pfb_f = NULL;
            int count = 0;
            int pfb_status = bcf_get_info_float(pfb_reader->readers[0].header, pfb_record, AF_key.c_str(), &pfb_f, &count);
            if (pfb_status < 0 || count == 0)
            {
                // std::cout << "Field " << AF_key << " not found, or count is 0" << std::endl;
                continue;
            }
            double pfb = (double) pfb_f[0];
            free(pfb_f);

            // Continue if the population frequency is outside the threshold
            if (pfb <= MIN_PFB || pfb >= MAX_PFB)
            {
                continue;
            }

            // Add the population frequency to the SNP data
            if (snp_pfb_map.find(pos) == snp_pfb_map.end())
            {
                snp_pfb_map[pos] = pfb;
            } else {
                // Keep the larger population frequency
                if (pfb > snp_pfb_map[pos])
                {
                    snp_pfb_map[pos] = pfb;
                }
            }

            // if (test_count < 10)
            // {
            //     std::cout << "Population frequency for " << chr << ":" << pos << " = " << pfb << std::endl;
            //     test_count++;
            // }
        }
            // std::cout << "Record: " << pfb_record->pos << std::endl;
            // std::cout << "QUAL: " << pfb_record->qual << std::endl;

            //     // Skip if not a SNP
            //     if (!bcf_is_snp(pfb_record))
            //     {
            //         std::cout << "Skipping non-SNP at " << chr << ":" << pos << std::endl;
            //         continue;
            //     }

            //     // Get the population frequency for the SNP
            //     float *pfb_f = NULL;
            //     int count = 0;
            //     int pfb_status = bcf_get_info_float(pfb_header, pfb_record, AF_key.c_str(), &pfb_f, &count);
            //     if (pfb_status < 0 || count == 0)
            //     {
            //         std::cout << "Field " << AF_key << " not found, or count is 0" << std::endl;
            //         continue;
            //     }
            //     double pfb = (double) pfb_f[0];
            //     free(pfb_f);

            //     // Continue if the population frequency is outside the threshold
            //     if (pfb <= MIN_PFB || pfb >= MAX_PFB)
            //     {
            //         continue;
            //     }

            //     // Add the population frequency to the SNP data
            //     // snp_info.insertSNPPopulationFrequency(chr_no_prefix, pos, pfb);
            //     if (snp_pfb_map.find(pos) == snp_pfb_map.end())
            //     {
            //         snp_pfb_map[pos] = pfb;
            //     } else {
            //         // Keep the larger population frequency
            //         if (pfb > snp_pfb_map[pos])
            //         {
            //             snp_pfb_map[pos] = pfb;
            //         }
            //     }
            //     if (print_count < 10)
            //     {
            //         std::cout << "Population frequency for " << chr << ":" << pos << " = " << pfb << std::endl;
            //         print_count++;
            //     }        }
    }
    if (pfb_reader->errnum)
    {
        std::cerr << "ERROR: " <<bcf_sr_strerror(pfb_reader->errnum) << std::endl;
    }

    // Clean up
    // bcf_destroy(pfb_record);
    bcf_sr_destroy(pfb_reader);

    // // Open the population frequency file
    // std::cout << "Opening population frequency file: " << pfb_filepath << std::endl;
    // htsFile *pfb_file = hts_open(pfb_filepath.c_str(), "r");
    // if (!pfb_file)
    // {
    //     throw std::runtime_error("ERROR: Could not open population frequency file: " + pfb_filepath);
    // }

    // // Enable multi-threading
    // std::cout << "Setting number of threads to " << thread_count << std::endl;
    // hts_set_threads(pfb_file, thread_count);

    // // Read the header
    // std::cout << "Reading header from population frequency file..." << std::endl;
    // bcf_hdr_t *pfb_header = bcf_hdr_read(pfb_file);
    // if (!pfb_header)
    // {
    //     bcf_close(pfb_file);
    //     throw std::runtime_error("ERROR: Could not read header from population frequency file: " + pfb_filepath);
    // }

    // // Load the index
    // hts_idx_t *pfb_index = bcf_index_load(pfb_filepath.c_str());
    // if (!pfb_index)
    // {
    //     bcf_hdr_destroy(pfb_header);
    //     bcf_close(pfb_file);
    //     throw std::runtime_error("ERROR: Could not load index for population frequency file: " + pfb_filepath);
    // }

    // // Construct the region string
    // std::string region_str = chr + ":" + std::to_string(start_pos) + "-" + std::to_string(end_pos);
    // hts_itr_t *pfb_iter = bcf_itr_querys(pfb_index, pfb_header, region_str.c_str());
    // if (!pfb_iter)
    // {
    //     // Try using the other chromosome notation
    //     std::string alt_region_str = "chr" + chr + ":" + std::to_string(start_pos) + "-" + std::to_string(end_pos);
    //     pfb_iter = bcf_itr_querys(pfb_index, pfb_header, alt_region_str.c_str());
    //     if (!pfb_iter)
    //     {
    //         hts_idx_destroy(pfb_index);
    //         bcf_hdr_destroy(pfb_header);
    //         bcf_close(pfb_file);
    //         throw std::runtime_error("ERROR: Could not create iterator for region: " + alt_region_str);
    //     } else {
    //         region_str = alt_region_str;
    //         std::cout << "Successfully created iterator for region: " << region_str << std::endl;
    //     }
    //     // hts_idx_destroy(pfb_index);
    //     // bcf_hdr_destroy(pfb_header);
    //     // bcf_close(pfb_file);
    //     // throw std::runtime_error("ERROR: Could not create iterator for region: " + region_str);
    // }

    // // Set up the record
    // std::cout << "Initializing BCF record..." << std::endl;
    // bcf1_t *pfb_record = bcf_init();
    // if (!pfb_record)
    // {
    //     bcf_hdr_destroy(pfb_header);
    //     bcf_close(pfb_file);
    //     throw std::runtime_error("ERROR: Could not initialize BCF record.");
    // }

    // // Read the population frequencies for the region
    // std::cout << "[TEST] Reading population frequencies for region " << region_str << " (AF_key = " << AF_key << ")..." << std::endl;
    // int print_count = 0;
    // int test_count = 0;
    // while (bcf_itr_next(pfb_file, pfb_iter, pfb_record) >= 0)
    // {
    //     test_count++;
    //     // Get the chromosome and position
    //     // std::cout << "Reading record..." << std::endl;
    //     uint32_t pos = pfb_record->pos + 1;  // 0-based to 1-based

    //     // Skip if not a SNP
    //     if (!bcf_is_snp(pfb_record))
    //     {
    //         std::cout << "Skipping non-SNP at " << chr << ":" << pos << std::endl;
    //         continue;
    //     }

    //     // Get the population frequency for the SNP
    //     float *pfb_f = NULL;
    //     int count = 0;
    //     int pfb_status = bcf_get_info_float(pfb_header, pfb_record, AF_key.c_str(), &pfb_f, &count);
    //     if (pfb_status < 0 || count == 0)
    //     {
    //         std::cout << "Field " << AF_key << " not found, or count is 0" << std::endl;
    //         continue;
    //     }
    //     double pfb = (double) pfb_f[0];
    //     free(pfb_f);

    //     // Continue if the population frequency is outside the threshold
    //     if (pfb <= MIN_PFB || pfb >= MAX_PFB)
    //     {
    //         continue;
    //     }

    //     // Add the population frequency to the SNP data
    //     // snp_info.insertSNPPopulationFrequency(chr_no_prefix, pos, pfb);
    //     if (snp_pfb_map.find(pos) == snp_pfb_map.end())
    //     {
    //         snp_pfb_map[pos] = pfb;
    //     } else {
    //         // Keep the larger population frequency
    //         if (pfb > snp_pfb_map[pos])
    //         {
    //             snp_pfb_map[pos] = pfb;
    //         }
    //     }
    //     if (print_count < 10)
    //     {
    //         std::cout << "Population frequency for " << chr << ":" << pos << " = " << pfb << std::endl;
    //         print_count++;
    //     }
    // }
    // std::cout << "Finished reading population frequencies for region " << region_str << std::endl;
    // std::cout << "Test count: " << test_count << std::endl;
}

void CNVCaller::saveSVCopyNumberToTSV(SNPData& snp_data, std::string filepath, std::string chr, uint32_t start, uint32_t end, std::string sv_type, double likelihood)
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
        uint32_t pos        = snp_data.pos[i];
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

void CNVCaller::updateSNPData(SNPData& snp_data, uint32_t pos, double pfb, double baf, double log2_cov, bool is_snp)
{
    // Update the SNP data
    snp_data.pos.emplace_back(pos);
    snp_data.pfb.emplace_back(pfb);
    snp_data.baf.emplace_back(baf);
    snp_data.log2_cov.emplace_back(log2_cov);
    snp_data.is_snp.emplace_back(is_snp);
}

void CNVCaller::querySNPs(std::string chr, uint32_t start, uint32_t end, std::set<uint32_t>& snp_pos, std::unordered_map<uint32_t, double>& snp_baf, std::unordered_map<uint32_t, double>& snp_pfb)
{
    std::string snp_chr = chr;
    chr = removeChrPrefix(chr);

    // Create an ordered map of SNP positions to BAF and PFB values
    std::map<uint32_t, std::tuple<double, double>> snp_map;

    // Query SNPs within a range (start, end) and return their BAF and PFB
    // values as separate vectors
    // std::vector<double> bafs;
    // std::vector<double> pfbs;
    // std::vector<uint32_t> pos;
    double pfb_default = 0.5;

    // Read the SNP data from the VCF file
    this->readSNPAlleleFrequencies(snp_chr, start, end, snp_pos, snp_baf);

    // Query the SNPs within the range and return their BAFs and corresponding
    // positions
    // auto snp_start = this->snp_baf_keys.lower_bound(start);
    // auto snp_end = this->snp_baf_keys.upper_bound(end);
    // if (snp_start == this->snp_baf_keys.end())
    // {
    //     // return std::make_tuple(pos, bafs, pfbs);
    //     return;
    // }

    // Query the population frequencies for the SNPs
    std::unordered_map<uint32_t, double> pfb_map;
    this->readSNPPopulationFrequencies(chr, start, end, pfb_map);

    // Filter out the SNP population frequencies that are not in the SNP
    // position set
    // std::unordered_map<uint32_t, double> snp_pfb;
    for (auto& pos : snp_pos)
    {
        if (pfb_map.find(pos) != pfb_map.end())
        {
            snp_pfb[pos] = pfb_map[pos];
        } else {
            snp_pfb[pos] = pfb_default;
        }
    }

    // // Get the PFB values for the SNPs from the keys
    // // Create the PFB vector using the SNP positions (loop through snp_pos,
    // // query the pfb_map, and push the value to the vector)
    // for (size_t i = 0; i < snp_pos.size(); i++)
    // {
    //     uint32_t snp_pos = snp_pos[i];
    //     double pfb = pfb_default;
    //     if (pfb_map.find(snp_pos) != pfb_map.end())
    //     {
    //         pfb = pfb_map[snp_pos];
    //     } else {
    //         pfb = pfb_default;
    //     }
    //     snp_pfb.push_back(pfb);
    // }

    // // Get the PFB values for the SNPs from the keys
    // for (auto it = snp_start; it != snp_end; it++)
    // {
    //     uint32_t snp_pos = *it;
    //     pos.push_back(snp_pos);
    //     bafs.push_back(this->snp_baf_map[snp_pos]);

    //     // Get the PFB value for the SNP
    //     if (this->snp_pfb_map.find(snp_pos) != this->snp_pfb_map.end())
    //     {
    //         pfbs.push_back(this->snp_pfb_map[snp_pos]);
    //     } else {
    //         pfbs.push_back(pfb_default);
    //     }
    // }
    // auto& baf_bst = this->snp_baf_map[chr];
    // auto baf_start = baf_bst.lower_bound({start, 0.0});
    // auto baf_end = baf_bst.upper_bound({end, 0.0});
    // for (auto it = baf_start; it != baf_end; it++) {
    //     bafs.push_back(std::get<1>(*it));
    //     pos.push_back(std::get<0>(*it));
    // }



    // auto& pfb_map = this->snp_pfb_map[chr];
    // for (size_t i = 0; i < pos.size(); i++) {
    //     uint32_t snp_pos = pos[i];
    //     if (pfb_map.find(snp_pos) != pfb_map.end()) {
    //         pfbs[i] = pfb_map[snp_pos];
    //     }
    // }
    
    // return std::make_tuple(pos, bafs, pfbs);
}
