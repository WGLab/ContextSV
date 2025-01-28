
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
    int region_length = (int) (end_pos - start_pos + 1);
    if (region_length < sample_size)
    {
        sample_size = region_length;
    }

    std::vector<uint32_t> snp_pos(sample_size, 0);
    std::vector<double> snp_baf(sample_size, -1.0);
    std::vector<double> snp_pfb(sample_size, 0.5);
    std::vector<double> snp_log2_cov(sample_size, 0.0);
    std::vector<bool> is_snp(sample_size, false);
    this->readSNPAlleleFrequencies(chr, start_pos, end_pos, snp_pos, snp_baf, snp_pfb, is_snp, input_data);

    // Get the log2 ratio for <sample_size> evenly spaced positions in the
    // region
    this->calculateRegionLog2Ratio(start_pos, end_pos, sample_size, pos_depth_map, mean_chr_cov, snp_log2_cov);

    // Update the SNP data with all information
    snp_data.pos = std::move(snp_pos);
    snp_data.baf = std::move(snp_baf);
    snp_data.pfb = std::move(snp_pfb);
    snp_data.log2_cov = std::move(snp_log2_cov);
    snp_data.is_snp = std::move(is_snp);
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
    if (input_data.getSaveCNVData())
    {
        uint32_t sv_half_length = (end_pos - start_pos) / 2.0;
        snp_start_pos = start_pos > sv_half_length ? start_pos - sv_half_length : 1;
        snp_end_pos = end_pos + sv_half_length;
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
    // double pct_threshold = 0.90;
    // double pct_threshold = 0.80;
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

    // Save the SV calls as a TSV file if enabled
    bool copy_number_change = (predicted_cnv_type != SVType::UNKNOWN && predicted_cnv_type != SVType::NEUTRAL);
    // if (save_cnv_data && copy_number_change && (end_pos - start_pos) > 10000)
    if (input_data.getSaveCNVData() && copy_number_change && (end_pos - start_pos) > 10000)
    {
        std::string cnv_type_str = getSVTypeString(predicted_cnv_type);
        const std::string output_dir = input_data.getOutputDir();
        std::string sv_filename = output_dir + "/" + cnv_type_str + "_" + chr + "_" + std::to_string((int) start_pos) + "-" + std::to_string((int) end_pos) + "_SPLITALN.tsv";
        printMessage("Saving SV split-alignment copy number predictions to " + sv_filename + "...");
        this->saveSVCopyNumberToTSV(snp_data, sv_filename, chr, start_pos, end_pos, cnv_type_str, likelihood);
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

        // Update the SV calls with the CNV type and genotype
        SVType updated_sv_type = getSVTypeFromCNState(max_state);
        std::string genotype = cnv_genotype_map.at(max_state);

        // Determine the SV calling method used to call the SV
        std::string data_type;
        data_type = "HMM";

        // Update the SV genotype if known
        // printMessage("Updating SV call for " + chr + ":" + std::to_string((int)start_pos) + "-" + std::to_string((int)end_pos) + "...");
        if (updated_sv_type != SVType::UNKNOWN)
        {
            sv_call.genotype = genotype;
            sv_call.data_type = data_type;
            sv_call.hmm_likelihood = likelihood;
        }

        // Update the SV type if known
        if (updated_sv_type != SVType::UNKNOWN && updated_sv_type != SVType::NEUTRAL)
        {
            std::string sv_type_str = getSVTypeString(updated_sv_type);
            sv_call.sv_type = sv_type_str;
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
double CNVCaller::calculateMeanChromosomeCoverage(std::string chr, std::vector<uint32_t>& chr_pos_depth_map, const std::string& bam_filepath, int thread_count) const
{
    {
        // Open the BAM file
        std::lock_guard<std::mutex> lock(this->shared_mutex);  // Lock the BAM file
        samFile *bam_file = sam_open(bam_filepath.c_str(), "r");
        if (!bam_file)
        {
            printError("ERROR: Could not open BAM file: " + bam_filepath);
            return 0.0;
        }

        // Enable multi-threading. This is possible here due to the lock
        hts_set_threads(bam_file, thread_count);
        // if (single_chr)
        // {
        //     hts_set_threads(bam_file, thread_count);
        // }

        // Read the header
        bam_hdr_t *bam_header = sam_hdr_read(bam_file);
        if (!bam_header)
        {
            sam_close(bam_file);
            printError("ERROR: Could not read header from BAM file: " + bam_filepath);
            return 0.0;
        }

        // Load the index
        hts_idx_t *bam_index = sam_index_load(bam_file, bam_filepath.c_str());
        if (!bam_index)
        {
            bam_hdr_destroy(bam_header);
            sam_close(bam_file);
            printError("ERROR: Could not load index for BAM file: " + bam_filepath);
            return 0.0;
        }
        BamFileGuard bam_guard(bam_file, bam_index, bam_header);  // Guard to close the BAM file

        // Create an iterator for the chromosome
        hts_itr_t *bam_iter = sam_itr_querys(bam_index, bam_header, chr.c_str());
        if (!bam_iter)
        {
            printError("ERROR: Could not create iterator for chromosome: " + chr + ", check if the chromosome exists in the BAM file.");
            return 0.0;
        }

        // Initialize the record
        bam1_t *bam_record = bam_init1();
        if (!bam_record)
        {
            hts_itr_destroy(bam_iter);
            printError("ERROR: Could not initialize BAM record.");
            return 0.0;
        }

        // Iterate through the chromosome and update the depth map
        while (sam_itr_next(bam_file, bam_iter, bam_record) >= 0)
        {
            // Ignore UNMAP, SECONDARY, QCFAIL, and DUP reads
            if (bam_record->core.flag & BAM_FUNMAP || bam_record->core.flag & BAM_FSECONDARY || bam_record->core.flag & BAM_FQCFAIL || bam_record->core.flag & BAM_FDUP)
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
                        try {
                            chr_pos_depth_map[ref_pos + j]++;
                        } catch (const std::out_of_range& oor) {
                            printError("Out of range error for " + chr + ":" + std::to_string(ref_pos+j));
                        }
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

        // Clean up
        bam_destroy1(bam_record);
        hts_itr_destroy(bam_iter);
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

    return mean_chr_cov;
}

void CNVCaller::calculateRegionLog2Ratio(uint32_t start_pos, uint32_t end_pos, int sample_size, const std::vector<uint32_t>& pos_depth_map, double mean_chr_cov, std::vector<double>& log2_region) const
{
    uint32_t region_length = end_pos - start_pos + 1;
    for (int i = 0; i < sample_size; i++)
    {
        uint32_t pos = start_pos + ((double)region_length / sample_size) * i;
        try {
            uint32_t depth = pos_depth_map.at(pos);

            // Calculate the log2 ratio for the position
            if (depth == 0)
            {
                log2_region[i] = 0.0;
            } else {
                log2_region[i] = log2((double) depth / mean_chr_cov);
            }

        } catch (const std::out_of_range& e) {
            log2_region[i] = 0.0;
        }
    }
}

void CNVCaller::readSNPAlleleFrequencies(std::string chr, uint32_t start_pos, uint32_t end_pos, std::vector<uint32_t>& snp_pos, std::vector<double>& snp_baf, std::vector<double>& snp_pfb, std::vector<bool>& is_snp, const InputData& input_data) const
{
    // Lock during reading
    std::lock_guard<std::mutex> lock(this->shared_mutex);

    // --------- SNP file ---------
    const std::string snp_filepath = input_data.getSNPFilepath();
    if (snp_filepath.empty())
    {
        printError("ERROR: SNP file path is empty.");
        return;
    }

    // Initialize the SNP file reader
    // printMessage("Initializing SNP reader...");
    bcf_srs_t *snp_reader = bcf_sr_init();
    if (!snp_reader)
    {
        printError("ERROR: Could not initialize SNP reader.");
        return;
    }
    snp_reader->require_index = 1;

    // Use multi-threading. This is possible here due to the lock
    int thread_count = input_data.getThreadCount();
    bcf_sr_set_threads(snp_reader, thread_count);

    // Add the SNP file to the reader
    // printMessage("Adding SNP file to reader...");
    if (bcf_sr_add_reader(snp_reader, snp_filepath.c_str()) < 0)
    {
        bcf_sr_destroy(snp_reader);
        printError("ERROR: Could not add SNP file to reader: " + snp_filepath);
        return;
    }
    // printMessage("SNP file added to reader.");

    // --------- Population allele frequency file ---------

    // Get the population allele frequency file path
    bool use_pfb = true;
    const std::string pfb_filepath = input_data.getAlleleFreqFilepath(chr);
    if (pfb_filepath.empty())
    {
        use_pfb = false;
    }
    
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
        // printMessage("Adding population allele frequency file to reader...");
        if (bcf_sr_add_reader(pfb_reader, pfb_filepath.c_str()) < 0)
        {
            printError("ERROR: Could not add population allele frequency file to reader: " + pfb_filepath);

            // Clean up
            bcf_sr_destroy(pfb_reader);
            bcf_sr_destroy(snp_reader);
            return;
        }

        // Use multi-threading. This is possible here due to the lock
        bcf_sr_set_threads(pfb_reader, thread_count);
    }

    // Split the region into samples
    int sample_size = snp_pos.size();
    std::vector<std::string> region_chunks = splitRegionIntoChunks(chr, start_pos, end_pos, sample_size);

    // Loop through the samples and read the SNP data, storing the first
    // SNP position and BAF value for each sample
    // int print_count = 0;
    int current_region = 0;
    for (size_t i = 0; i < region_chunks.size(); ++i)
    {
        current_region++;
        // Lock during reading
        // std::lock_guard<std::mutex> lock(this->shared_mutex);

        // Read the SNP data ----------------------------------------------

        // Set the region
        // printMessage("Setting region for SNP reader...");
        std::string region_str = region_chunks[i];
        if (bcf_sr_set_regions(snp_reader, region_str.c_str(), 0) < 0)
        {
            printError("ERROR: Could not set region for SNP reader: " + region_str);
            break;
        }
        // printMessage("Region set for SNP reader, loading SNP data...");

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
                snp_pos[i] = pos;
                snp_baf[i] = baf;
                is_snp[i] = true;
                snp_found = true;

                break;  // Only one SNP per region
            }
        }

        if (snp_reader->errnum)
        {
            printError("ERROR: " + std::string(bcf_sr_strerror(snp_reader->errnum)));
        }

        // Continue if no SNP was found in the region
        if (!snp_found)
        {
            continue;
        }

        // Read the population allele frequency data ----------------------
        if (use_pfb)
        {
            // Set the region as the SNP position
            // printMessage("Setting region for population allele frequency reader...");
            uint32_t target_snp_pos = snp_pos[i];  // Already 1-based
            std::string snp_region_str = chr_gnomad + ":" + std::to_string(target_snp_pos) + "-" + std::to_string(target_snp_pos);
            if (bcf_sr_set_regions(pfb_reader, snp_region_str.c_str(), 0) < 0)
            {
                printError("ERROR: Could not set region for population allele frequency reader: " + region_str);
                break;
            }
            // printMessage("Region set for population allele frequency reader, loading population allele frequency data...");

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

                // if (!bcf_sr_has_line(pfb_reader, 0))
                // {
                //     continue;
                // }
                // bcf1_t *pfb_record = bcf_sr_get_line(pfb_reader, 0);
                // if (pfb_record)
                // {
                //     // Skip if not a SNP
                //     if (!bcf_is_snp(pfb_record))
                //     {
                //         continue;
                //     }

                // Get the population frequency for the SNP
                // float *pfb_f = NULL;
                // int count = 0;
                int pfb_status = bcf_get_info_float(pfb_reader->readers[0].header, pfb_record, AF_key.c_str(), &pfb_f, &count);
                if (pfb_status < 0 || count == 0)
                {
                    continue;
                }
                // double pfb = (double) pfb_f[0];
                double pfb = static_cast<double>(pfb_f[0]);
                // free(pfb_f);

                // Skip if outside the acceptable range
                if (pfb <= MIN_PFB || pfb >= MAX_PFB)
                {
                    continue;
                }

                // Add the population frequency to the SNP data
                snp_pfb[i] = pfb;

                // Break after finding the SNP position
                break;

                // if (print_count < 20) {
                //     printMessage("SNP " + std::to_string(snp_pos[i]) + " BAF: " + std::to_string(snp_baf[i]) + " PFB: " + std::to_string(snp_pfb[i]) + " (Region: " + snp_region_str + ")");
                //     print_count++;
                // }
            }
            free(pfb_f);

            // }
            if (pfb_reader->errnum)
            {
                printError("ERROR: " + std::string(bcf_sr_strerror(pfb_reader->errnum)));
            }
        }
        // printMessage("SNP region " + std::to_string(current_region) + " of " + std::to_string(region_chunks.size()) + " completed.");
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

void CNVCaller::updateSNPData(SNPData& snp_data, uint32_t pos, double pfb, double baf, double log2_cov, bool is_snp)
{
    // Update the SNP data
    snp_data.pos.emplace_back(pos);
    snp_data.pfb.emplace_back(pfb);
    snp_data.baf.emplace_back(baf);
    snp_data.log2_cov.emplace_back(log2_cov);
    snp_data.is_snp.emplace_back(is_snp);
}
