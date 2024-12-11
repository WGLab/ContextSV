
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
#include <utility>    // std::pair

#include "utils.h"
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
void CNVCaller::runViterbi(const CHMM& hmm, SNPData& snp_data, std::pair<std::vector<int>, double>& prediction)
{
    int data_count = (int) snp_data.pos.size();
    if (data_count == 0)
    {
        // throw std::runtime_error("Error: No SNP data found for Viterbi
        // algorithm.");
        printError("ERROR: No SNP data found for Viterbi algorithm.");
        prediction = std::make_pair(std::vector<int>(), 0.0);
    }
    prediction = testVit_CHMM(hmm, data_count, snp_data.log2_cov, snp_data.baf, snp_data.pfb);
}

// Function to obtain SNP information for a region
std::pair<SNPData, bool> CNVCaller::querySNPRegion(std::string chr, uint32_t start_pos, uint32_t end_pos, const std::vector<uint32_t>& pos_depth_map, double mean_chr_cov)
{
    SNPData snp_data;
    bool snps_found = false;
    uint32_t window_size = (uint32_t)this->input_data.getWindowSize();

    // Query the SNPs for the entire region
    std::set<uint32_t> snp_pos;
    std::unordered_map<uint32_t, double> snp_baf;
    std::unordered_map<uint32_t, double> snp_pfb;
    this->querySNPs(chr, start_pos, end_pos, snp_pos, snp_baf, snp_pfb);

    // Loop through the range of the SV region and query the SNPs in a sliding
    // window, then calculate the log2 ratio for each window
    for (uint32_t i = start_pos; i <= end_pos; i += window_size)
    {
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

                // Calculate the log2 ratio for the SNP bin
                double bin_cov = calculateLog2Ratio(bin_start, bin_end, pos_depth_map, mean_chr_cov);
                this->updateSNPData(snp_data, snp_window_pos[j], snp_window_pfbs[j], snp_window_bafs[j], bin_cov, true);

                // Update the previous bin start
                bin_start = bin_end + 1;
            }
        }
    }

    return std::make_pair(snp_data, snps_found);
}

std::tuple<double, SVType, std::string, bool> CNVCaller::runCopyNumberPrediction(std::string chr, const CHMM& hmm, uint32_t start_pos, uint32_t end_pos, double mean_chr_cov, const std::vector<uint32_t>& pos_depth_map)
{
    // Check that the start position is less than the end position
    if (start_pos >= end_pos)
    {
        // throw std::runtime_error("ERROR: Invalid SV region for copy number prediction: " + chr + ":" + std::to_string((int)start_pos) + "-" + std::to_string((int)end_pos));
        printError("ERROR: Invalid SV region for copy number prediction: " + chr + ":" + std::to_string((int)start_pos) + "-" + std::to_string((int)end_pos));
        return std::make_tuple(0.0, SVType::UNKNOWN, "./.", false);
    }

    // Run the Viterbi algorithm on SNPs in the SV region +/- 1/2
    // the SV length
    // Only extened the region if "save CNV data" is enabled
    uint32_t snp_start_pos = start_pos;
    uint32_t snp_end_pos = end_pos;
    if (this->input_data.getSaveCNVData())
    {
        uint32_t sv_half_length = (end_pos - start_pos) / 2.0;
        snp_start_pos = start_pos > sv_half_length ? start_pos - sv_half_length : 1;
        snp_end_pos = end_pos + sv_half_length;
    }
    // uint32_t sv_half_length = (end_pos - start_pos) / 2.0;
    // uint32_t snp_start_pos = start_pos > sv_half_length ? start_pos - sv_half_length : 1;
    // uint32_t snp_end_pos = end_pos + sv_half_length;

    // Query the SNP region for the SV candidate
    std::pair<SNPData, bool> snp_call = querySNPRegion(chr, snp_start_pos, snp_end_pos, pos_depth_map, mean_chr_cov);
    SNPData& sv_snps = snp_call.first;
    bool sv_snps_found = snp_call.second;

    // Run the Viterbi algorithm
    // printMessage("Running Viterbi algorithm for SV " + chr + ":" + std::to_string((int)start_pos) + "-" + std::to_string((int)end_pos) + " (" + std::to_string(sv_snps.pos.size()) + " SNPs, start=" + std::to_string(snp_start_pos) + ", end=" + std::to_string(snp_end_pos) + ")...");
    std::pair<std::vector<int>, double> prediction;
    runViterbi(hmm, sv_snps, prediction);
    if (prediction.first.size() == 0)
    {
        return std::make_tuple(0.0, SVType::UNKNOWN, "./.", sv_snps_found);
    }

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
        printMessage("Saving SV split-alignment copy number predictions to " + sv_filename + "...");
        this->saveSVCopyNumberToTSV(sv_snps, sv_filename, chr, start_pos, end_pos, cnv_type_str, likelihood);
    }
    
    return std::make_tuple(likelihood, predicted_cnv_type, genotype, sv_snps_found);
}


void CNVCaller::runCIGARCopyNumberPrediction(std::string chr, std::vector<SVCall> &sv_candidates, const CHMM& hmm, double mean_chr_cov, const std::vector<uint32_t>& pos_depth_map)
{
    // Map with counts for each CNV type
    std::map<int, int> cnv_type_counts;
    for (int i = 0; i < 6; i++)
    {
        cnv_type_counts[i] = 0;
    }
    
    // Loop through each SV candidate and predict the copy number state
    int min_length = this->input_data.getMinCNVLength();
    for (auto& sv_call : sv_candidates)
    {

        // Get the SV candidate
        uint32_t start_pos = sv_call.start;
        uint32_t end_pos = sv_call.end;
        
        // Error if start > end
        if (start_pos >= end_pos)
        {
        	std::cerr << "Position error for CIGAR SV at " << chr << ":" << start_pos << "-" << end_pos << std::endl;
        	continue;
        }

        // Skip if not the minimum length for CNV predictions
        if ((end_pos - start_pos) < (uint32_t) min_length)
        {
            continue;
        }

        // Loop through the SV region +/- 1/2 SV length and run copy number
        // predictions
        // Only extend the region if "save CNV data" is enabled
        uint32_t snp_start_pos = start_pos;
        uint32_t snp_end_pos = end_pos;
        if (this->input_data.getSaveCNVData())
        {
            uint32_t sv_half_length = (end_pos - start_pos) / 2.0;
            snp_start_pos = start_pos > sv_half_length ? start_pos - sv_half_length : 1;
            snp_end_pos = end_pos + sv_half_length;
        }
        std::pair<SNPData, bool> snp_call = this->querySNPRegion(chr, snp_start_pos, snp_end_pos, pos_depth_map, mean_chr_cov);
        SNPData& sv_snps = snp_call.first;
        bool snps_found = snp_call.second;

        // Run the Viterbi algorithm
        if (sv_snps.pos.size() == 0) {
        	std::cerr << "ERROR: No windows for SV " << chr << ":" << start_pos << "-" << end_pos << " (" << snp_start_pos << "," << snp_end_pos << std::endl;
        	continue;
        }
        
        std::pair<std::vector<int>, double> prediction;
        runViterbi(hmm, sv_snps, prediction);
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

        // Update the SV genotype if known
        if (updated_sv_type != SVType::UNKNOWN)
        {
            sv_call.genotype = genotype;
            sv_call.data_type = data_type;
            sv_call.hmm_likelihood = likelihood;
        }

        // Update the SV type if known
        // printMessage("Updating SV copy number data for SV " + chr + ":" + std::to_string((int)start_pos) + "-" + std::to_string((int)end_pos) + "...");
        if (updated_sv_type != SVType::UNKNOWN && updated_sv_type != SVType::NEUTRAL)
        {
            std::string sv_type_str = getSVTypeString(updated_sv_type);
            sv_call.sv_type = sv_type_str;
            // sv_call.data_type = data_type;
            // sv_call.genotype = genotype;
            // sv_call.hmm_likelihood = likelihood;
        }

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
            printMessage("Saving SV CIGAR copy number predictions to " + sv_filename);
            this->saveSVCopyNumberToTSV(sv_snps, sv_filename, chr, start_pos, end_pos, cnv_type_str, likelihood);
        }
    }
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

// Calculate the mean chromosome coverage
std::pair<double, std::vector<uint32_t>> CNVCaller::calculateMeanChromosomeCoverage(std::string chr, uint32_t chr_len)
{
    std::vector<uint32_t> chr_pos_depth_map(chr_len+1, 0); // 1-based index
    {
        // Lock the bam file
        std::lock_guard<std::mutex> lock(this->bam_file_mtx);

        // Open the BAM file
        std::string bam_filepath = this->input_data.getShortReadBam();
        samFile *bam_file = sam_open(bam_filepath.c_str(), "r");
        if (!bam_file)
        {
            // throw std::runtime_error("ERROR: Could not open BAM file: " +
            // bam_filepath);
            printError("ERROR: Could not open BAM file: " + bam_filepath);
            return std::make_pair(0.0, chr_pos_depth_map);
        }

        // Enable multi-threading
        // hts_set_threads(bam_file, this->input_data.getThreadCount());

        // Read the header
        bam_hdr_t *bam_header = sam_hdr_read(bam_file);
        if (!bam_header)
        {
            sam_close(bam_file);
            printError("ERROR: Could not read header from BAM file: " + bam_filepath);
            return std::make_pair(0.0, chr_pos_depth_map);
            // throw std::runtime_error("ERROR: Could not read header from BAM file: " + bam_filepath);
        }

        // Load the index
        hts_idx_t *bam_index = sam_index_load(bam_file, bam_filepath.c_str());
        if (!bam_index)
        {
            bam_hdr_destroy(bam_header);
            sam_close(bam_file);
            // throw std::runtime_error("ERROR: Could not load index for BAM
            // file: " + bam_filepath);
            printError("ERROR: Could not load index for BAM file: " + bam_filepath);
            return std::make_pair(0.0, chr_pos_depth_map);  
        }

        // Create an iterator for the chromosome
        hts_itr_t *bam_iter = sam_itr_querys(bam_index, bam_header, chr.c_str());
        if (!bam_iter)
        {
            hts_idx_destroy(bam_index);
            bam_hdr_destroy(bam_header);
            sam_close(bam_file);
            // throw std::runtime_error("ERROR: Could not create iterator for
            // chromosome: " + chr + ", check if the chromosome exists in the
            // BAM file.");
            printError("ERROR: Could not create iterator for chromosome: " + chr + ", check if the chromosome exists in the BAM file.");
            return std::make_pair(0.0, chr_pos_depth_map);
        }

        // Initialize the record
        bam1_t *bam_record = bam_init1();
        if (!bam_record)
        {
            hts_itr_destroy(bam_iter);
            hts_idx_destroy(bam_index);
            bam_hdr_destroy(bam_header);
            sam_close(bam_file);
            // throw std::runtime_error("ERROR: Could not initialize BAM
            // record.");
            printError("ERROR: Could not initialize BAM record.");
            return std::make_pair(0.0, chr_pos_depth_map);
        }

        // Iterate through the chromosome and update the depth map
        // std::unordered_map<uint32_t, int> chr_pos_depth_map;
        while (sam_itr_next(bam_file, bam_iter, bam_record) >= 0)
        {
            // Ignore UNMAP, SECONDARY, QCFAIL, and DUP reads
            if (bam_record->core.flag & BAM_FUNMAP || bam_record->core.flag & BAM_FSECONDARY || bam_record->core.flag & BAM_FQCFAIL || bam_record->core.flag & BAM_FDUP)
            {
                continue;
            }
            
            // Parse the CIGAR string to get the depth (match, sequence match, and
            // mismatch)
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
                        try {
                            chr_pos_depth_map[ref_pos + j]++;
                        } catch (const std::out_of_range& oor) {
                            // std::cerr << "Out of range error for " << chr <<
                            // ":" << ref_pos+j << std::endl;
                            printError("Out of range error for " + chr + ":" + std::to_string(ref_pos+j));
                        }
                        // chr_pos_depth_map[ref_pos + j]++;
                    }
                }
                
                // Update the reference coordinate based on the CIGAR operation
                // https://samtools.github.io/hts-specs/SAMv1.pdf
                if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP || op == BAM_CEQUAL || op == BAM_CDIFF) {
                    ref_pos += op_len;
                } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP || op == BAM_CPAD) {
                    // Do nothing
                } else {
                    // throw std::runtime_error("ERROR: Unknown CIGAR operation:
                    // " + std::to_string(op));
                    printError("ERROR: Unknown CIGAR operation: " + std::to_string(op));
                }
            }
        }

        // Clean up
        bam_destroy1(bam_record);
        hts_itr_destroy(bam_iter);
        hts_idx_destroy(bam_index);
        bam_hdr_destroy(bam_header);
        sam_close(bam_file);
    }

    // Calculate the mean chromosome coverage for positions with non-zero depth
    uint64_t cum_depth = 0;
    uint32_t pos_count = 0;
    for (const auto& pos_depth : chr_pos_depth_map)
    {
        if (pos_depth > 0)
        {
            cum_depth += pos_depth;
            pos_count++;
        }
    }

    double mean_chr_cov = 0.0;
    if (pos_count > 0)
    {
        mean_chr_cov = static_cast<double>(cum_depth) / static_cast<double>(pos_count);
    }

    return std::make_pair(mean_chr_cov, chr_pos_depth_map);
}

double CNVCaller::calculateLog2Ratio(uint32_t start_pos, uint32_t end_pos, const std::vector<uint32_t>& pos_depth_map, double mean_chr_cov)
{
    // Use the position and depth map to calculate the log2 ratio
    double cum_depth = 0;
    int pos_count = 0;
    for (uint32_t i = start_pos; i <= end_pos; i++)
    {
        if (i < pos_depth_map.size() && pos_depth_map[i] > 0)
        {
            cum_depth += pos_depth_map[i];
            pos_count++;
        }
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
        // throw std::runtime_error("ERROR: SNP file path is empty.");
        printError("ERROR: SNP file path is empty.");
        return;
    }

    // Initialize the synced reader
    bcf_srs_t *snp_reader = bcf_sr_init();
    if (!snp_reader)
    {
        // throw std::runtime_error("ERROR: Could not initialize SNP reader.");
        printError("ERROR: Could not initialize SNP reader.");
        return;
    }

    // Lock during reading
    std::lock_guard<std::mutex> lock(this->snp_file_mtx);

    // Set the region
    std::string region_str = chr + ":" + std::to_string(start_pos) + "-" + std::to_string(end_pos);
    if (bcf_sr_set_regions(snp_reader, region_str.c_str(), 0) < 0)
    {
        bcf_sr_destroy(snp_reader);
        // throw std::runtime_error("ERROR: Could not set region for SNP reader:
        // " + region_str);
        printError("ERROR: Could not set region for SNP reader: " + region_str);
        return;
    }

    // Set multi-threading
    // int thread_count = this->input_data.getThreadCount();
    // bcf_sr_set_threads(snp_reader, thread_count);

    // Enable index usage
    snp_reader->require_index = 1;

    // Add the SNP file to the reader
    if (bcf_sr_add_reader(snp_reader, snp_filepath.c_str()) < 0)
    {
        bcf_sr_destroy(snp_reader);
        // throw std::runtime_error("ERROR: Could not add SNP file to reader: "
        // + snp_filepath);
        printError("ERROR: Could not add SNP file to reader: " + snp_filepath);
        return;
    }

    // Get the header
    bcf_hdr_t *snp_header = bcf_sr_get_header(snp_reader, 0);
    if (!snp_header)
    {
        bcf_sr_destroy(snp_reader);
        // throw std::runtime_error("ERROR: Could not get header for SNP
        // reader.");
        printError("ERROR: Could not get header for SNP reader.");
        return;
    }

    // std::cout << "Iterating through SNPs in region " << region_str << "..." << std::endl;
    int record_count = 0;
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
                // throw std::runtime_error("ERROR: AD value is missing for SNP
                // at " + chr + ":" + std::to_string(pos));
                printError("ERROR: AD value is missing for SNP at " + chr + ":" + std::to_string(pos));
                continue;
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
        }
    }

    // Clean up
    bcf_sr_destroy(snp_reader);
}

void CNVCaller::readSNPPopulationFrequencies(std::string chr, uint32_t start_pos, uint32_t end_pos, std::unordered_map<uint32_t, double>& snp_pfb_map)
{
    // Get the population frequency file for the chromosome
    std::string pfb_filepath = this->input_data.getAlleleFreqFilepath(chr);
    if (pfb_filepath == "")
    {
        // printError("No population frequency file provided for chromosome " + chr);
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
    // int thread_count = this->input_data.getThreadCount();

    // Initialize the synced reader
    bcf_srs_t *pfb_reader = bcf_sr_init();
    if (!pfb_reader)
    {
        // throw std::runtime_error("ERROR: Could not initialize synced reader
        // for population frequency file: " + pfb_filepath);
        printError("ERROR: Could not initialize synced reader for population frequency file: " + pfb_filepath);
        return;
    }

    // Lock during reading
    std::lock_guard<std::mutex> lock(this->pfb_file_mtx);

    // Set the region for the synced reader
    std::string region_str = chr_gnomad + ":" + std::to_string(start_pos) + "-" + std::to_string(end_pos);
    if (bcf_sr_set_regions(pfb_reader, region_str.c_str(), 0) < 0)
    {
        bcf_sr_destroy(pfb_reader);
        // throw std::runtime_error("ERROR: Could not set region for synced
        // reader: " + region_str);
        printError("ERROR: Could not set region for synced reader: " + region_str);
        return;
    }

    // Set multi-threading
    // bcf_sr_set_threads(pfb_reader, thread_count);

    // Enable index usage
    pfb_reader->require_index = 1;

    // Add the population frequency file to the synced reader
    if (bcf_sr_add_reader(pfb_reader, pfb_filepath.c_str()) < 0)
    {
        bcf_sr_destroy(pfb_reader);
        // throw std::runtime_error("ERROR: Could not add population frequency
        // file to synced reader: " + pfb_filepath);
        printError("ERROR: Could not add population frequency file to synced reader: " + pfb_filepath);
        return;
    }

    // Get the header
    bcf_hdr_t *pfb_header = bcf_sr_get_header(pfb_reader, 0);
    if (!pfb_header)
    {
        bcf_sr_destroy(pfb_reader);
        // throw std::runtime_error("ERROR: Could not get header for population
        // frequency file: " + pfb_filepath);
        printError("ERROR: Could not get header for population frequency file: " + pfb_filepath);
        return;
    }

    int record_count = 0;
    while (bcf_sr_next_line(pfb_reader) > 0)
    {
        if (!bcf_sr_has_line(pfb_reader, 0))
        {
            continue;
        }
        // pfb_record = bcf_sr_get_line(pfb_reader, 0);
        bcf1_t *pfb_record = bcf_sr_get_line(pfb_reader, 0);
        // Do something with the record
        if (pfb_record)
        {
            record_count++;
            // Skip if not a SNP
            if (!bcf_is_snp(pfb_record))
            {
                continue;
            }

            uint32_t pos = (uint32_t) pfb_record->pos + 1;  // 0-based to 1-based

            // Get the population frequency for the SNP
            float *pfb_f = NULL;
            int count = 0;
            int pfb_status = bcf_get_info_float(pfb_reader->readers[0].header, pfb_record, AF_key.c_str(), &pfb_f, &count);
            if (pfb_status < 0 || count == 0)
            {
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
        }
    }
    if (pfb_reader->errnum)
    {
        // std::cerr << "ERROR: " <<bcf_sr_strerror(pfb_reader->errnum) <<
        // std::endl;
        printError("ERROR: " + std::string(bcf_sr_strerror(pfb_reader->errnum)));
    }

    // Clean up
    bcf_sr_destroy(pfb_reader);
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

    // Query the SNP allele frequencies for the SNPs
    std::map<uint32_t, std::tuple<double, double>> snp_map;
    this->readSNPAlleleFrequencies(snp_chr, start, end, snp_pos, snp_baf);

    // Query the population frequencies for the SNPs
    std::unordered_map<uint32_t, double> pfb_map;
    this->readSNPPopulationFrequencies(chr, start, end, pfb_map);

    // Filter out the SNP population frequencies that are not in the SNP
    // position set
    double pfb_default = 0.5;
    for (auto& pos : snp_pos)
    {
        if (pfb_map.find(pos) != pfb_map.end())
        {
            snp_pfb[pos] = pfb_map[pos];
        } else {
            snp_pfb[pos] = pfb_default;
        }
    }
}
