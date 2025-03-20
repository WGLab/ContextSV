
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
#include <unordered_set>
#include <execution>  // std::execution::par

#include "utils.h"
#include "sv_types.h"

#define MIN_PFB 0.01
#define MAX_PFB 0.99
/// @endcond

using namespace sv_types;


// Function to call the Viterbi algorithm for the CHMM
void CNVCaller::runViterbi(const CHMM& hmm, SNPData& snp_data, std::pair<std::vector<int>, double>& prediction) const
{
    int data_count = (int) snp_data.pos.size();
    if (data_count == 0)
    {
        printError("ERROR: No SNP data found for Viterbi algorithm.");
        prediction = std::make_pair(std::vector<int>(), 0.0);
    }
    prediction = testVit_CHMM(hmm, data_count, snp_data.log2_cov, snp_data.baf, snp_data.pfb);
}

// Function to obtain SNP information for a region
void CNVCaller::querySNPRegion(std::string chr, uint32_t start_pos, uint32_t end_pos, const std::vector<uint32_t>& pos_depth_map, double mean_chr_cov, SNPData& snp_data, const InputData& input_data) const
{
    // Initialize the SNP data with default values and sample size length
    int sample_size = input_data.getSampleSize();
    std::vector<uint32_t> snp_pos;
    std::unordered_map<uint32_t, double> snp_baf_map;
    std::unordered_map<uint32_t, double> snp_pfb_map;
    this->readSNPAlleleFrequencies(chr, start_pos, end_pos, snp_pos, snp_baf_map, snp_pfb_map, input_data);

    // Get the log2 ratio for <sample_size> evenly spaced positions in the
    // region
    sample_size = std::max((int) snp_pos.size(), sample_size);

    // Loop through evenly spaced positions in the region and get the log2 ratio
    double pos_step = (double) (end_pos - start_pos + 1) / (double) sample_size;
    std::unordered_set<uint32_t> snp_pos_set(snp_pos.begin(), snp_pos.end());
    std::unordered_map<std::string, double> window_log2_map;
    for (int i = 0; i < sample_size; i++)
    {
        uint32_t window_start = (uint32_t) (start_pos + i * pos_step);
        uint32_t window_end = (uint32_t) (start_pos + (i + 1) * pos_step);

        // Calculate the mean depth for the window
        double cov_sum = 0.0;
        int pos_count = 0;
        for (int j = 0; j < pos_step; j++)
        {
            uint32_t pos = (uint32_t) (start_pos + i * pos_step + j);
            if (pos > end_pos)
            {
                break;
            }
            if (pos < pos_depth_map.size()) {
                cov_sum += pos_depth_map[pos];
                pos_count++;
            }

        }
        double log2_cov = 0.0;
        if (pos_count > 0)
        {
            if (cov_sum == 0)
            {
                // Use a small value to avoid division by zero
                cov_sum = 1e-9;
            }
            log2_cov = log2((cov_sum / (double) pos_count) / mean_chr_cov);
        }

        // Store the log2 ratio for the window
        std::string window_key = std::to_string(window_start) + "-" + std::to_string(window_end);
        window_log2_map[window_key] = log2_cov;
    }

    // Create new vectors for the SNP data
    std::vector<uint32_t> snp_pos_hmm;
    std::vector<double> snp_baf_hmm;
    std::vector<double> snp_pfb_hmm;
    std::vector<double> snp_log2_hmm;
    std::vector<bool> is_snp_hmm;

    // Loop through the window ranges and append all SNPs in the range, using
    // the log2 ratio for the window
    for (const auto& window : window_log2_map)
    {
        uint32_t window_start = std::stoi(window.first.substr(0, window.first.find('-')));
        uint32_t window_end = std::stoi(window.first.substr(window.first.find('-') + 1));
        double log2_cov = window.second;

        // Loop through the SNP positions and add them to the SNP data
        bool snp_found = false;
        for (uint32_t pos : snp_pos)
        {
            if (pos >= window_start && pos <= window_end)
            {
                snp_pos_hmm.push_back(pos);
                snp_baf_hmm.push_back(snp_baf_map[pos]);
                snp_pfb_hmm.push_back(snp_pfb_map[pos]);
                snp_log2_hmm.push_back(log2_cov);
                is_snp_hmm.push_back(true);
                snp_found = true;
            }
        }
        if (!snp_found)
        {
            // If no SNPs were found in the window, add a dummy SNP with the
            // log2 ratio for the window, using the window center as the SNP
            // position
            uint32_t window_center = (window_start + window_end) / 2;
            snp_pos_hmm.push_back(window_center);
            snp_baf_hmm.push_back(-1.0);
            snp_pfb_hmm.push_back(0.5);
            snp_log2_hmm.push_back(log2_cov);
            is_snp_hmm.push_back(false);
        }
    }

    // Update the SNP data with all information
    snp_data.pos = std::move(snp_pos_hmm);
    snp_data.baf = std::move(snp_baf_hmm);
    snp_data.pfb = std::move(snp_pfb_hmm);
    snp_data.log2_cov = std::move(snp_log2_hmm);
    snp_data.is_snp = std::move(is_snp_hmm);
}

std::tuple<double, SVType, std::string, bool> CNVCaller::runCopyNumberPrediction(std::string chr, const CHMM& hmm, uint32_t start_pos, uint32_t end_pos, double mean_chr_cov, const std::vector<uint32_t>& pos_depth_map, const InputData& input_data) const
{
    // Check that the start position is less than the end position
    if (start_pos > end_pos)
    {
        printError("ERROR: Invalid SV region for copy number prediction: " + chr + ":" + std::to_string((int)start_pos) + "-" + std::to_string((int)end_pos));
        return std::make_tuple(0.0, SVType::UNKNOWN, "./.", false);
    }

    // Run the Viterbi algorithm on SNPs in the SV region
    // Only extend the region if "save CNV data" is enabled
    uint32_t snp_start_pos = start_pos;
    uint32_t snp_end_pos = end_pos;
    SNPData before_sv;
    SNPData after_sv;
    if (input_data.getSaveCNVData())
    {
        uint32_t sv_half_length = (end_pos - start_pos) / 2.0;
        if (start_pos > 1)
        {
            uint32_t before_sv_start = std::max((uint32_t) 1, start_pos - sv_half_length);
            uint32_t before_sv_end = start_pos - 1;
            querySNPRegion(chr, before_sv_start, before_sv_end, pos_depth_map, mean_chr_cov, before_sv, input_data);
        }
        uint32_t chr_last_index = pos_depth_map.size() - 1;
        if (end_pos < chr_last_index)
        {
            uint32_t after_sv_start = end_pos + 1;
            uint32_t after_sv_end = std::min(chr_last_index, end_pos + sv_half_length);
            querySNPRegion(chr, after_sv_start, after_sv_end, pos_depth_map, mean_chr_cov, after_sv, input_data);
        }
    }

    // Query the SNP region for the SV candidate
    SNPData snp_data;
    querySNPRegion(chr, snp_start_pos, snp_end_pos, pos_depth_map, mean_chr_cov, snp_data, input_data);

    // Run the Viterbi algorithm
    std::pair<std::vector<int>, double> prediction;
    runViterbi(hmm, snp_data, prediction);
    if (prediction.first.size() == 0)
    {
        return std::make_tuple(0.0, SVType::UNKNOWN, "./.", false);
    }

    std::vector<int>& state_sequence = prediction.first;
    double likelihood = prediction.second;

    // Get all the states in the SV region
    std::vector<int> sv_states;
    for (size_t i = 0; i < state_sequence.size(); i++)
    {
        if (snp_data.pos[i] >= start_pos && snp_data.pos[i] <= end_pos)
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
        genotype = cnv_genotype_map.at(max_state);
    }
    snp_data.state_sequence = std::move(state_sequence);  // Move the state sequence to the SNP data

    // Save the SV calls if enabled
    bool copy_number_change = (predicted_cnv_type != SVType::UNKNOWN && predicted_cnv_type != SVType::NEUTRAL);
    if (input_data.getSaveCNVData() && copy_number_change && (end_pos - start_pos) > 50000)
    {
        // Set B-allele and population frequency values to 0 for non-SNPs
        for (size_t i = 0; i < snp_data.pos.size(); i++)
        {
            if (!snp_data.is_snp[i])
            {
                snp_data.baf[i] = 0.0;
                snp_data.pfb[i] = 0.0;
            }
        }
        for (size_t i = 0; i < before_sv.pos.size(); i++)
        {
            if (!before_sv.is_snp[i])
            {
                before_sv.baf[i] = 0.0;
                before_sv.pfb[i] = 0.0;
            }
        }
        for (size_t i = 0; i < after_sv.pos.size(); i++)
        {
            if (!after_sv.is_snp[i])
            {
                after_sv.baf[i] = 0.0;
                after_sv.pfb[i] = 0.0;
            }
        }

        // Save the SNP data to JSON
        std::string cnv_type_str = getSVTypeString(predicted_cnv_type);
        std::string json_filepath = input_data.getCNVOutputFile();
        printMessage("Saving SV copy number predictions to " + json_filepath + "...");

        this->saveSVCopyNumberToJSON(before_sv, after_sv, snp_data, chr, start_pos, end_pos, cnv_type_str, likelihood, json_filepath);
    }
    
    return std::make_tuple(likelihood, predicted_cnv_type, genotype, true);
}


void CNVCaller::runCIGARCopyNumberPrediction(std::string chr, std::vector<SVCall>& sv_candidates, const CHMM& hmm, double mean_chr_cov, const std::vector<uint32_t>& pos_depth_map, const InputData& input_data) const
{
    // Map with counts for each CNV type
    std::map<int, int> cnv_type_counts;
    for (int i = 0; i < 6; i++)
    {
        cnv_type_counts[i] = 0;
    }
    
    // Loop through each SV candidate and predict the copy number state
    for (auto& sv_call : sv_candidates)
    {

        // Get the SV candidate
        uint32_t start_pos = sv_call.start;
        uint32_t end_pos = sv_call.end;
        
        // Error if start > end
        if (start_pos > end_pos)
        {
            printError("ERROR: Invalid SV region for copy number prediction: " + chr + ":" + std::to_string((int)start_pos) + "-" + std::to_string((int)end_pos));
        	continue;
        }

        // Skip if not the minimum length for CNV predictions
        if ((end_pos - start_pos) < input_data.getMinCNVLength())
        {
            continue;
        }

        // Only extend the region if "save CNV data" is enabled
        SNPData snp_data;
        this->querySNPRegion(chr, start_pos, end_pos, pos_depth_map, mean_chr_cov, snp_data, input_data);

        // Run the Viterbi algorithm
        if (snp_data.pos.size() == 0) {
            printError("ERROR: No SNP data found for Viterbi algorithm for CIGAR SV at " + chr + ":" + std::to_string((int)start_pos) + "-" + std::to_string((int)end_pos));
        	continue;
        }
        
        std::pair<std::vector<int>, double> prediction;
        runViterbi(hmm, snp_data, prediction);
        std::vector<int>& state_sequence = prediction.first;
        double likelihood = prediction.second;

        // Get all the states in the SV region
        std::vector<int> sv_states;
        for (size_t i = 0; i < state_sequence.size(); i++)
        {
            if (snp_data.pos[i] >= start_pos && snp_data.pos[i] <= end_pos)
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

        // Update the SV information if it does not conflict with the current SV type
        SVType updated_sv_type = getSVTypeFromCNState(max_state);
        bool is_valid_update = isValidCopyNumberUpdate(sv_call.sv_type, updated_sv_type);
        if (is_valid_update)
        {
            std::string genotype = cnv_genotype_map.at(max_state);
            std::string data_type = "CIGAR+HMM";
            sv_call.sv_type = updated_sv_type;
            sv_call.hmm_likelihood = likelihood;
            sv_call.genotype = genotype;
            sv_call.data_type = data_type;
        }
    }
}

std::vector<std::string> CNVCaller::splitRegionIntoChunks(std::string chr, uint32_t start_pos, uint32_t end_pos, int chunk_count) const
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
void CNVCaller::calculateMeanChromosomeCoverage(const std::vector<std::string>& chromosomes, std::unordered_map<std::string, std::vector<uint32_t>>& chr_pos_depth_map, std::unordered_map<std::string, double>& chr_mean_cov_map, const std::string& bam_filepath, int thread_count) const
{
    // Open the BAM file
    // std::shared_lock<std::shared_mutex> lock(this->shared_mutex);  // Lock the BAM file
    printMessage("Opening BAM file: " + bam_filepath);
    samFile *bam_file = sam_open(bam_filepath.c_str(), "r");
    if (!bam_file)
    {
        printError("ERROR: Could not open BAM file: " + bam_filepath);
        return;
    }

    // Enable multi-threading while opening the BAM file
    hts_set_threads(bam_file, thread_count);

    // Read the header
    bam_hdr_t *bam_header = sam_hdr_read(bam_file);
    if (!bam_header)
    {
        sam_close(bam_file);
        printError("ERROR: Could not read header from BAM file: " + bam_filepath);
        return;
    }

    // Load the index
    hts_idx_t *bam_index = sam_index_load(bam_file, bam_filepath.c_str());
    if (!bam_index)
    {
        bam_hdr_destroy(bam_header);
        sam_close(bam_file);
        printError("ERROR: Could not load index for BAM file: " + bam_filepath);
        return;
    }
    BamFileGuard bam_guard(bam_file, bam_index, bam_header);  // Guard to close the BAM file

    // Initialize the record
    bam1_t *bam_record = bam_init1();
    if (!bam_record)
    {
        // Clean up the BAM file and index
        bam_hdr_destroy(bam_header);
        sam_close(bam_file);
        printError("ERROR: Could not initialize BAM record.");
        return;
    }

    // Iterate through each chromosome and update the depth map
    int current_chr = 0;
    int total_chr_count = chromosomes.size();
    for (const std::string& chr : chromosomes)
    {
        // Create an iterator for the chromosome
        hts_itr_t *bam_iter = sam_itr_querys(bam_index, bam_header, chr.c_str());
        if (!bam_iter)
        {
            printError("ERROR: Could not create iterator for chromosome: " + chr + ", check if the chromosome exists in the BAM file.");
            continue;
        }

        printMessage("(" + std::to_string(++current_chr) + "/" + std::to_string(total_chr_count) + ") Reading BAM file for chromosome: " + chr);
        std::vector<uint32_t>& pos_depth_map = chr_pos_depth_map[chr];
        while (sam_itr_next(bam_file, bam_iter, bam_record) >= 0)
        {
            // Ignore UNMAP, SECONDARY, QCFAIL, and DUP reads
            uint16_t flag = bam_record->core.flag;
            if (flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP))
            {
                continue;
            }

            // Parse the CIGAR string to get the depth (match, sequence match, and
            // mismatch)
            uint32_t pos = (uint32_t)bam_record->core.pos + 1;  // 0-based to 1-based
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
                        if (ref_pos + j >= pos_depth_map.size())
                        {
                            printError("ERROR: Reference position out of range for " + chr + ":" + std::to_string(ref_pos+j));
                            continue;
                        }
                        pos_depth_map[ref_pos + j]++;
                    }
                }
                
                // Update the reference coordinate based on the CIGAR operation
                // https://samtools.github.io/hts-specs/SAMv1.pdf
                if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP || op == BAM_CEQUAL || op == BAM_CDIFF) {
                    ref_pos += op_len;
                } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP || op == BAM_CPAD) {
                    // Do nothing
                } else {
                    printError("ERROR: Unknown CIGAR operation: " + std::to_string(op));
                }
            }
        }
        hts_itr_destroy(bam_iter);
        
        // Parallel sum of the depth map
        uint64_t cum_depth = std::reduce(
            std::execution::par,
            pos_depth_map.begin(),
            pos_depth_map.end(),
            0ULL
        );

        // Parallel count of the non-zero depth positions
        uint32_t pos_count = std::count_if(
            std::execution::par,
            pos_depth_map.begin(),
            pos_depth_map.end(),
            [](uint32_t depth) { return depth > 0; }
        );

        double mean_chr_cov = (pos_count > 0) ? static_cast<double>(cum_depth) / static_cast<double>(pos_count) : 0.0;
        chr_mean_cov_map[chr] = mean_chr_cov;
    }
}

void CNVCaller::readSNPAlleleFrequencies(std::string chr, uint32_t start_pos, uint32_t end_pos, std::vector<uint32_t>& snp_pos, std::unordered_map<uint32_t, double>& snp_baf, std::unordered_map<uint32_t, double>& snp_pfb, const InputData& input_data) const
{
    // Lock during reading
    std::shared_lock<std::shared_mutex> lock(this->shared_mutex);

    // --------- SNP file ---------
    const std::string snp_filepath = input_data.getSNPFilepath();
    if (snp_filepath.empty())
    {
        printError("ERROR: SNP file path is empty.");
        return;
    }

    // Initialize the SNP file reader
    bcf_srs_t *snp_reader = bcf_sr_init();
    if (!snp_reader)
    {
        printError("ERROR: Could not initialize SNP reader.");
        return;
    }
    snp_reader->require_index = 1;

    // Use multi-threading if not threading by chromosome
    int thread_count = input_data.getThreadCount();
    bcf_sr_set_threads(snp_reader, thread_count);

    // Add the SNP file to the reader
    if (bcf_sr_add_reader(snp_reader, snp_filepath.c_str()) < 0)
    {
        bcf_sr_destroy(snp_reader);
        printError("ERROR: Could not add SNP file to reader: " + snp_filepath);
        return;
    }

    // --------- Population allele frequency file ---------

    // Get the population allele frequency file path
    bool use_pfb = true;
    const std::string pfb_filepath = input_data.getAlleleFreqFilepath(chr);
    if (pfb_filepath.empty())
    {
        use_pfb = false;
    }

    // Ensure the file exists (ifsstream will throw an exception if the file
    // does not exist)
    std::ifstream pfb_file(pfb_filepath);
    if (!pfb_file)
    {
        use_pfb = false;
    }
    pfb_file.close();

    bcf_srs_t *pfb_reader = bcf_sr_init();
    std::string chr_gnomad;
    std::string AF_key;
    if (use_pfb)
    {
        // Determine the ethnicity-specific allele frequency key
        AF_key = "AF";
        const std::string eth = input_data.getEthnicity();
        if (eth != "")
        {
            AF_key += "_" + eth;
        }

        // Check if the filepath uses the 'chr' prefix notations based on the
        // chromosome name (*.chr1.vcf.gz vs *.1.vcf.gz)
        chr_gnomad = chr;  // gnomAD data may or may not have the 'chr' prefix
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

        // Initialize the population allele frequency reader
        if (!pfb_reader)
        {
            printError("ERROR: Could not initialize population allele frequency reader.");

            // Clean up
            bcf_sr_destroy(snp_reader);
            return;
        }
        pfb_reader->require_index = 1;

        // Add the population allele frequency file to the reader
        if (bcf_sr_add_reader(pfb_reader, pfb_filepath.c_str()) < 0)
        {
            printError("ERROR: Could not add population allele frequency file to reader: " + pfb_filepath);

            // Clean up
            bcf_sr_destroy(pfb_reader);
            bcf_sr_destroy(snp_reader);
            return;
        }

        // Use multi-threading if not threading by chromosome
        int thread_count = input_data.getThreadCount();
        bcf_sr_set_threads(pfb_reader, thread_count);
    }

    // Read the SNP data ----------------------------------------------
    // Set the region
    std::string region_str = chr + ":" + std::to_string(start_pos) + "-" + std::to_string(end_pos);
    if (bcf_sr_set_regions(snp_reader, region_str.c_str(), 0) < 0)  //chr.c_str(), 0) < 0)
    {
        printError("ERROR: Could not set region for SNP reader: " + chr);
        bcf_sr_destroy(snp_reader);
        bcf_sr_destroy(pfb_reader);
        return;
    }

    bool snp_found = false;
    while (bcf_sr_next_line(snp_reader) > 0)
    {
        if (!bcf_sr_has_line(snp_reader, 0))
        {
            continue;
        }
        bcf1_t *snp_record = bcf_sr_get_line(snp_reader, 0);
        if (snp_record)
        {
            uint32_t pos = (uint32_t)snp_record->pos + 1;

            // Skip if not a SNP
            if (!bcf_is_snp(snp_record))
            {
                continue;
            }

            // Get the QUAL, DP, and AD values
            if (bcf_float_is_missing(snp_record->qual) || snp_record->qual <= 30)
            {
                continue;
            }

            // Extract DP from FORMAT field
            int32_t *dp = 0;
            int dp_count = 0;
            int dp_ret = bcf_get_format_int32(snp_reader->readers[0].header, snp_record, "DP", &dp, &dp_count);
            if (dp_ret < 0 || dp[0] <= 10)
            {
                continue;
            }
            free(dp);

            // Skip if the SNP does not pass the filter
            if (bcf_has_filter(snp_reader->readers[0].header, snp_record, const_cast<char*>("PASS")) != 1)
            {
                continue;
            }

            // Extract AD from FORMAT field
            int32_t *ad = 0;
            int ad_count = 0;
            int ad_ret = bcf_get_format_int32(snp_reader->readers[0].header, snp_record, "AD", &ad, &ad_count);
            if (ad_ret < 0 || ad_count < 2)
            {
                continue;
            }

            // Calculate the B-allele frequency (BAF)
            double baf = (double) ad[1] / (double) (ad[0] + ad[1]);
            free(ad);

            // Add the SNP position and BAF information
            snp_pos.push_back(pos);
            snp_baf[pos] = baf;
            printMessage("SNP found: " + chr + ":" + std::to_string(pos) + " BAF: " + std::to_string(baf));
            snp_found = true;
        }
    }

    if (snp_reader->errnum)
    {
        printError("ERROR: " + std::string(bcf_sr_strerror(snp_reader->errnum)));
    }

    // Continue if no SNP was found in the region
    if (!snp_found)
    {
        bcf_sr_destroy(snp_reader);
        bcf_sr_destroy(pfb_reader);
        return;
    }

    // Read the population allele frequency data ----------------------
    // Get the minimum and maximum SNP positions
    uint32_t min_snp_pos = *std::min_element(snp_pos.begin(), snp_pos.end());
    uint32_t max_snp_pos = *std::max_element(snp_pos.begin(), snp_pos.end());
    std::unordered_set<uint32_t> snp_pos_set(snp_pos.begin(), snp_pos.end());
    if (use_pfb)
    {
        // Set the region for the population allele frequency reader
        std::string pfb_region_str = chr_gnomad + ":" + std::to_string(min_snp_pos) + "-" + std::to_string(max_snp_pos);
        if (bcf_sr_set_regions(pfb_reader, pfb_region_str.c_str(), 0) < 0)
        {
            printError("ERROR: Could not set region for population allele frequency reader: " + pfb_region_str);
        }

        // Find the SNP position in the population allele frequency file
        float *pfb_f = NULL;
        int count = 0;
        while (bcf_sr_next_line(pfb_reader) > 0)
        {
            // Get the SNP record and validate
            bcf1_t *pfb_record = bcf_sr_get_line(pfb_reader, 0);
            if (!pfb_record || !bcf_is_snp(pfb_record))
            {
                continue;  // Skip if not a SNP
            }

            // Get the SNP position
            uint32_t pfb_pos = (uint32_t) pfb_record->pos + 1;
            if (snp_pos_set.find(pfb_pos) == snp_pos_set.end())
            {
                continue;  // Skip if the SNP position is not in the set
            }

            // Get the population frequency for the SNP
            int pfb_status = bcf_get_info_float(pfb_reader->readers[0].header, pfb_record, AF_key.c_str(), &pfb_f, &count);
            if (pfb_status < 0 || count == 0)
            {
                continue;
            }
            double pfb = static_cast<double>(pfb_f[0]);

            // Skip if outside the acceptable range
            if (pfb <= MIN_PFB || pfb >= MAX_PFB)
            {
                continue;
            }
            // snp_pfb[i] = pfb;
            snp_pfb[pfb_pos] = pfb;
            printMessage("Population frequency found: " + chr + ":" + std::to_string(pfb_pos) + " PFB: " + std::to_string(pfb));
            break;
        }
        free(pfb_f);
    }
    
    // Clean up
    bcf_sr_destroy(snp_reader);
    bcf_sr_destroy(pfb_reader);
}

void CNVCaller::saveSVCopyNumberToTSV(SNPData& snp_data, std::string filepath, std::string chr, uint32_t start, uint32_t end, std::string sv_type, double likelihood) const
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

void CNVCaller::saveSVCopyNumberToJSON(SNPData &before_sv, SNPData &after_sv, SNPData &snp_data, std::string chr, uint32_t start, uint32_t end, std::string sv_type, double likelihood, const std::string& filepath) const
{
    // Append the SV information to the JSON file
    std::ofstream json_file(filepath, std::ios::app);
    if (!json_file.is_open())
    {
        std::cerr << "ERROR: Could not open JSON file for writing: " << filepath << std::endl;
        exit(1);
    }

    // If not the first record, write the closing bracket
    // Check if file is empty
    if (isFileEmpty(filepath))
    {
        json_file << "[\n";
    } else {
        // Close the previous JSON object
        json_file << "},\n";
    }

    json_file << "{\n";
    json_file << "  \"chromosome\": \"" << chr << "\",\n";
    json_file << "  \"start\": " << start << ",\n";
    json_file << "  \"end\": " << end << ",\n";
    json_file << "  \"sv_type\": \"" << sv_type << "\",\n";
    json_file << "  \"likelihood\": " << likelihood << ",\n";
    json_file << "  \"size\": " << (end - start + 1) << ",\n";
    json_file << "  \"before_sv\": {\n";
    json_file << "    \"positions\": [";
        for (size_t i = 0; i < before_sv.pos.size(); ++i)
        {
            json_file << before_sv.pos[i];
            if (i < before_sv.pos.size() - 1)
                json_file << ", ";
        }
        json_file << "],\n";
    json_file << "    \"b_allele_freq\": [";
        for (size_t i = 0; i < before_sv.baf.size(); ++i)
        {
            json_file << before_sv.baf[i];
            if (i < before_sv.baf.size() - 1)
                json_file << ", ";
        }
        json_file << "],\n";
    json_file << "    \"population_freq\": [";
        for (size_t i = 0; i < before_sv.pfb.size(); ++i)
        {
            json_file << before_sv.pfb[i];
            if (i < before_sv.pfb.size() - 1)
                json_file << ", ";
        }
        json_file << "],\n";
    json_file << "    \"log2_ratio\": [";
        for (size_t i = 0; i < before_sv.log2_cov.size(); ++i)
        {
            json_file << before_sv.log2_cov[i];
            if (i < before_sv.log2_cov.size() - 1)
                json_file << ", ";
        }
        json_file << "],\n";
    json_file << "    \"is_snp\": [";
        for (size_t i = 0; i < before_sv.is_snp.size(); ++i)
        {
            json_file << before_sv.is_snp[i];
            if (i < before_sv.is_snp.size() - 1)
                json_file << ", ";
        }
        json_file << "]\n";
    json_file << "  },\n";
    json_file << "  \"after_sv\": {\n";
    json_file << "    \"positions\": [";
        for (size_t i = 0; i < after_sv.pos.size(); ++i)
        {
            json_file << after_sv.pos[i];
            if (i < after_sv.pos.size() - 1)
                json_file << ", ";
        }
        json_file << "],\n";
    json_file << "    \"b_allele_freq\": [";
        for (size_t i = 0; i < after_sv.baf.size(); ++i)
        {
            json_file << after_sv.baf[i];
            if (i < after_sv.baf.size() - 1)
                json_file << ", ";
        }
        json_file << "],\n";
    json_file << "    \"population_freq\": [";
        for (size_t i = 0; i < after_sv.pfb.size(); ++i)
        {
            json_file << after_sv.pfb[i];
            if (i < after_sv.pfb.size() - 1)
                json_file << ", ";
        }
        json_file << "],\n";
    json_file << "    \"log2_ratio\": [";
        for (size_t i = 0; i < after_sv.log2_cov.size(); ++i)
        {
            json_file << after_sv.log2_cov[i];
            if (i < after_sv.log2_cov.size() - 1)
                json_file << ", ";
        }
        json_file << "],\n";
    json_file << "    \"is_snp\": [";
        for (size_t i = 0; i < after_sv.is_snp.size(); ++i)
        {
            json_file << after_sv.is_snp[i];
            if (i < after_sv.is_snp.size() - 1)
                json_file << ", ";
        }
        json_file << "]\n";
    json_file << "  },\n";
    json_file << "  \"sv\": {\n";
    json_file << "    \"positions\": [";
        for (size_t i = 0; i < snp_data.pos.size(); ++i)
        {
            json_file << snp_data.pos[i];
            if (i < snp_data.pos.size() - 1)
                json_file << ", ";
        }
        json_file << "],\n";
    json_file << "    \"b_allele_freq\": [";
        for (size_t i = 0; i < snp_data.baf.size(); ++i)
        {
            json_file << snp_data.baf[i];
            if (i < snp_data.baf.size() - 1)
                json_file << ", ";
        }
        json_file << "],\n";
    json_file << "    \"population_freq\": [";
        for (size_t i = 0; i < snp_data.pfb.size(); ++i)
        {
            json_file << snp_data.pfb[i];
            if (i < snp_data.pfb.size() - 1)
                json_file << ", ";
        }
        json_file << "],\n";
    json_file << "    \"log2_ratio\": [";
        for (size_t i = 0; i < snp_data.log2_cov.size(); ++i)
        {
            json_file << snp_data.log2_cov[i];
            if (i < snp_data.log2_cov.size() - 1)
                json_file << ", ";
        }
        json_file << "],\n";
    json_file << "    \"states\": [";
        for (size_t i = 0; i < snp_data.state_sequence.size(); ++i)
        {
            json_file << snp_data.state_sequence[i];
            if (i < snp_data.state_sequence.size() - 1)
                json_file << ", ";
        }
        json_file << "],\n";
    json_file << "    \"is_snp\": [";
        for (size_t i = 0; i < snp_data.is_snp.size(); ++i)
        {
            json_file << snp_data.is_snp[i];
            if (i < snp_data.is_snp.size() - 1)
                json_file << ", ";
        }
        json_file << "]\n";
    json_file << "  }\n";
    // json_file << "},\n";
    json_file.close();
    printMessage("Saved copy number predictions for " + chr + ":" + std::to_string(start) + "-" + std::to_string(end) + " to " + filepath);
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
