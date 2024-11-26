//
// sv_caller.cpp:
// Detect SVs from long read alignments
//

#include "sv_caller.h"

#include <htslib/sam.h>

/// @cond
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <thread>
#include <chrono>
#include <future>
#include <cmath>
#include <algorithm>
#include <fstream>

#include "utils.h"
#include "sv_types.h"
/// @endcond

# define DUP_SEQSIM_THRESHOLD 0.9  // Sequence similarity threshold for duplication detection

SVCaller::SVCaller(InputData &input_data)
    : input_data(input_data)  // Initialize the input data
{
}

int SVCaller::readNextAlignment(samFile *fp_in, hts_itr_t *itr, bam1_t *bam1)
{
    int ret = sam_itr_next(fp_in, itr, bam1);
    return ret;
}

// RegionData SVCaller::detectSVsFromRegion(std::string region)
std::tuple<std::set<SVCall>, PrimaryMap, SuppMap> SVCaller::detectCIGARSVs(samFile* fp_in, hts_idx_t* idx, bam_hdr_t* bamHdr, const std::string& region)
{
    // Create a read and iterator for the region
    bam1_t *bam1 = bam_init1();
    if (!bam1) {
        hts_idx_destroy(idx);
        bam_hdr_destroy(bamHdr);
        sam_close(fp_in);
        throw std::runtime_error("ERROR: failed to initialize BAM record");
    }
    hts_itr_t *itr = sam_itr_querys(idx, bamHdr, region.c_str());
    if (!itr) {
        bam_destroy1(bam1);
        hts_idx_destroy(idx);
        bam_hdr_destroy(bamHdr);
        sam_close(fp_in);
        throw std::runtime_error("ERROR: failed to query region " + region);
    }

    // Main loop to process the alignments
    std::set<SVCall> sv_calls;
    int num_alignments = 0;
    PrimaryMap primary_alignments;
    SuppMap supplementary_alignments;
    while (readNextAlignment(fp_in, itr, bam1) >= 0) {

        // Skip secondary and unmapped alignments, duplicates, QC failures, and low mapping quality
        if (bam1->core.flag & BAM_FSECONDARY || bam1->core.flag & BAM_FUNMAP || bam1->core.flag & BAM_FDUP || bam1->core.flag & BAM_FQCFAIL || bam1->core.qual < this->min_mapq) {
            continue;
        }
        const std::string qname = bam_get_qname(bam1);  // Query template name

        // Process primary alignments
        if (!(bam1->core.flag & BAM_FSUPPLEMENTARY)) {

            // Get the primary alignment information
            std::string chr = bamHdr->target_name[bam1->core.tid];
            int64_t start = bam1->core.pos;
            int64_t end = bam_endpos(bam1);  // This is the first position after the alignment
            bool fwd_strand = !(bam1->core.flag & BAM_FREVERSE);

            // Call SVs directly from the CIGAR string
            std::tuple<std::unordered_map<int, int>, int32_t, int32_t> query_info = this->detectSVsFromCIGAR(bamHdr, bam1, sv_calls, true);
            std::unordered_map<int, int> match_map = std::get<0>(query_info);
            int32_t query_start = std::get<1>(query_info);
            int32_t query_end = std::get<2>(query_info);

            // Add the primary alignment to the map
            AlignmentData alignment(chr, start, end, ".", query_start, query_end, match_map, fwd_strand);
            primary_alignments[qname] = alignment;

        // Process supplementary alignments
        } else if (bam1->core.flag & BAM_FSUPPLEMENTARY) {

            // Get the supplementary alignment information
            std::string chr = bamHdr->target_name[bam1->core.tid];
            int32_t start = bam1->core.pos;
            int32_t end = bam_endpos(bam1);
            bool fwd_strand = !(bam1->core.flag & BAM_FREVERSE);

            // Get CIGAR string information, but don't call SVs
            std::tuple<std::unordered_map<int, int>, int32_t, int32_t> query_info = this->detectSVsFromCIGAR(bamHdr, bam1, sv_calls, false);
            const std::unordered_map<int, int>& match_map = std::get<0>(query_info);
            int32_t query_start = std::get<1>(query_info);
            int32_t query_end = std::get<2>(query_info);

            // Add the supplementary alignment to the map
            AlignmentData alignment(chr, start, end, ".", query_start, query_end, match_map, fwd_strand);
            supplementary_alignments[qname].emplace_back(alignment);
        }

        num_alignments++;
    }

    // Clean up the iterator and alignment
    hts_itr_destroy(itr);
    bam_destroy1(bam1);

    return std::make_tuple(sv_calls, primary_alignments, supplementary_alignments);
}

double SVCaller::calculateMismatchRate(std::unordered_map<int, int> &match_map, int32_t start, int32_t end)
{
    int match_count = 0;
    int mismatch_count = 0;
    for (int i = start; i <= end; i++) {
        if (match_map.find(i) != match_map.end()) {
            if (match_map[i] == 1) {
                match_count++;
            } else {
                mismatch_count++;
            }
        }
    }
    double mismatch_rate = (double)mismatch_count / (double)(match_count + mismatch_count);

    return mismatch_rate;
}

std::tuple<std::unordered_map<int, int>, int32_t, int32_t> SVCaller::detectSVsFromCIGAR(bam_hdr_t* header, bam1_t* alignment, std::set<SVCall>& sv_calls, bool is_primary)
{
    std::string chr = header->target_name[alignment->core.tid];  // Chromosome name
    int32_t pos = alignment->core.pos;  // Leftmost position of the alignment in the reference genome (0-based)
    uint32_t* cigar = bam_get_cigar(alignment);  // CIGAR array
    int cigar_len = alignment->core.n_cigar;
    int query_pos = 0;
    std::unordered_map<int, int> query_match_map;  // Query position to match/mismatch (1/0) map

    // Loop through the CIGAR string, process operations, detect SVs (primary
    // only), update clipped base support, calculate sequence identity for
    // potential duplications (primary only), and calculate
    // the clipped base support and mismatch rate
    int32_t ref_pos;
    int32_t ref_end;
    int32_t query_start = 0;  // First alignment position in the query
    int32_t query_end = 0;    // Last alignment position in the query
    bool first_op = false;  // First alignment operation for the query
    double default_lh = 0.0;
    for (int i = 0; i < cigar_len; i++) {

        int op = bam_cigar_op(cigar[i]);  // CIGAR operation
        int op_len = bam_cigar_oplen(cigar[i]);  // CIGAR operation length
        
        // Process the CIGAR operation
        if (op == BAM_CINS && is_primary) {
            if (op_len >= this->min_sv_size) {

                // Get the sequence of the insertion from the query
                std::string ins_seq_str(op_len, ' ');
                for (int j = 0; j < op_len; j++) {
                    ins_seq_str[j] = seq_nt16_str[bam_seqi(bam_get_seq(alignment), query_pos + j)];
                }

                // To determine whether the insertion is a duplication, check
                // for sequence identity between the insertion and the
                // reference genome (duplications are typically >= 90%)

                // Loop through the reference sequence and calculate the
                // sequence identity +/- insertion length from the insertion
                // position.
                bool is_duplication = false;
                int ins_ref_pos;
                for (int j = pos - op_len; j <= pos; j++) {

                    // Get the string for the window (1-based coordinates)
                    ins_ref_pos = j + 1;
                    std::string window_str = this->input_data.queryRefGenome(chr, ins_ref_pos, ins_ref_pos + op_len - 1);

                    // Continue if the window string is empty (out-of-range)
                    if (window_str == "") {
                        continue;
                    }

                    // Calculate the sequence identity
                    int num_matches = 0;
                    for (int k = 0; k < op_len; k++) {
                        if (ins_seq_str[k] == window_str[k]) {
                            num_matches++;
                        }
                    }
                    float seq_identity = (float)num_matches / (float)op_len;

                    // Check if the target sequence identity is reached
                    if (seq_identity >= DUP_SEQSIM_THRESHOLD) {
                        is_duplication = true;
                        break;
                    }
                }

                // Determine whether to use a symbolic allele (>50bp) or the
                // actual sequence
                if (op_len > 50) {
                    ins_seq_str = "<INS>";
                } else {
                    ins_seq_str = ins_seq_str;
                }

                // Add to SV calls (1-based) with the appropriate SV type
                ref_pos = pos+1;
                ref_end = ref_pos + op_len -1;
                if (is_duplication) {
                    addSVCall(sv_calls, (uint32_t)ref_pos, (uint32_t)ref_end, "DUP", ins_seq_str, "CIGARDUP", "./.", default_lh);
                } else {
                    addSVCall(sv_calls, (uint32_t)ref_pos, (uint32_t)ref_end, "INS", ins_seq_str, "CIGARINS", "./.", default_lh);
                }
            }

        // Check if the CIGAR operation is a deletion
        } else if (op == BAM_CDEL && is_primary) {

            // Add the SV if greater than the minimum SV size
            if (op_len >= this->min_sv_size)
            {
                ref_pos = pos+1;
                ref_end = ref_pos + op_len -1;
                addSVCall(sv_calls, (uint32_t)ref_pos, (uint32_t)ref_end, "DEL", ".", "CIGARDEL", "./.", default_lh);
            }

        // Check if the CIGAR operation is a clipped base
        } else if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) {

            // sv_calls.updateClippedBaseSupport(chr, pos);  // Update clipped base support

            // Update the query alignment start position
            if (!first_op) {
                query_start = query_pos + op_len;
                first_op = true;
            }
        }

        // Update match/mismatch query map
        if (op == BAM_CEQUAL) {
            for (int j = 0; j < op_len; j++) {
                query_match_map[query_pos + j] = 1;
            }
        } else if (op == BAM_CDIFF) {
            for (int j = 0; j < op_len; j++) {
                query_match_map[query_pos + j] = 0;
            }
        } else if (op == BAM_CMATCH) {
            // Get the read sequence
            uint8_t* seq_ptr = bam_get_seq(alignment);
            std::string cmatch_seq_str = "";
            for (int j = 0; j < op_len; j++) {
                cmatch_seq_str += seq_nt16_str[bam_seqi(seq_ptr, query_pos + j)];
            }

            // Get the corresponding reference sequence
            int cmatch_pos = pos + 1;  // Querying the reference genome is 1-based
            std::string cmatch_ref_str = this->input_data.queryRefGenome(chr, cmatch_pos, cmatch_pos + op_len - 1);

            // Check that the two sequence lengths are equal
            if (cmatch_seq_str.length() != cmatch_ref_str.length()) {
                throw std::runtime_error("ERROR: Sequence lengths do not match for CIGAR operation: " + std::to_string(op));
            }

            // Compare the two sequences and update the mismatch map
            for (int j = 0; j < op_len; j++) {
                if (cmatch_seq_str[j] != cmatch_ref_str[j]) {
                    query_match_map[query_pos + j] = 0;
                } else {
                    query_match_map[query_pos + j] = 1;
                }
            }
        }

        // Update the reference coordinate based on the CIGAR operation
        // https://samtools.github.io/hts-specs/SAMv1.pdf
        if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP || op == BAM_CEQUAL || op == BAM_CDIFF) {
            pos += op_len;
        } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP || op == BAM_CPAD) {
            // Do nothing
        } else {
            throw std::runtime_error("ERROR: Unknown CIGAR operation: " + std::to_string(op));
        }

        // Update the query position based on the CIGAR operation (M, I, S, H)
        if (op == BAM_CMATCH || op == BAM_CINS || op == BAM_CSOFT_CLIP || op == BAM_CEQUAL || op == BAM_CDIFF) {
            query_pos += op_len;
        } else if (op == BAM_CDEL || op == BAM_CREF_SKIP || op == BAM_CHARD_CLIP || op == BAM_CPAD) {
            // Do nothing
        } else {
            throw std::runtime_error("ERROR: Unknown CIGAR operation: " + std::to_string(op));
        }
    }

    query_end = query_pos;  // Last alignment position in the query

    return std::tuple<std::unordered_map<int, int>, int32_t, int32_t>(query_match_map, query_start, query_end);
}

void SVCaller::processChromosome(const std::string& chr, const std::string& bam_filepath, CHMM hmm, std::set<SVCall>& combined_sv_calls, int min_cnv_length)
{
    // Open the BAM file
    samFile *fp_in = sam_open(bam_filepath.c_str(), "r");
    if (!fp_in) {
        throw std::runtime_error("ERROR: failed to open " + bam_filepath);
    }

    // Load the header
    bam_hdr_t *bamHdr = sam_hdr_read(fp_in);
    if (!bamHdr) {
        sam_close(fp_in);
        throw std::runtime_error("ERROR: failed to read header from " + bam_filepath);
    }

    // Load the index
    hts_idx_t *idx = sam_index_load(fp_in, bam_filepath.c_str());
    if (!idx) {
        bam_hdr_destroy(bamHdr);
        sam_close(fp_in);
        throw std::runtime_error("ERROR: failed to load index for " + bam_filepath);
    }

    // Split the chromosome into chunks for memory efficiency
    std::vector<std::string> region_chunks;
    int chunk_count = 100;
    if (this->input_data.isRegionSet()) {

        // Use one chunk for the specified region
        std::pair<int32_t, int32_t> region = this->input_data.getRegion();
        int region_start = region.first;
        int region_end = region.second;
        std::string chunk = chr + ":" + std::to_string(region_start) + "-" + std::to_string(region_end);
        region_chunks.push_back(chunk);
        // std::cout << "Using specified region " << chunk << "..." << std::endl;
        
    } else {
        int chr_len = this->input_data.getRefGenomeChromosomeLength(chr);
        int chunk_size = std::ceil((double)chr_len / chunk_count);
        for (int i = 0; i < chunk_count; i++) {
            int start = i * chunk_size + 1;  // 1-based
            int end = start + chunk_size;
            if (i == chunk_count - 1) {
                end = chr_len;
            }
            std::string chunk = chr + ":" + std::to_string(start) + "-" + std::to_string(end);
            region_chunks.push_back(chunk);
        }
        printMessage("Split chromosome " + chr + " into " + std::to_string(region_chunks.size()) + " chunks of size " + std::to_string(chunk_size) + "...");
    }

    // Load chromosome data for copy number predictions
    // std::cout << "Loading chromosome data for copy number predictions..." << std::endl;
    CNVCaller cnv_caller(this->input_data);
    cnv_caller.loadChromosomeData(chr);

    // Process each chunk one at a time
    // std::cout << "Processing " << region_chunks.size() << " region(s) for chromosome " << chr << "..." << std::endl;
    int region_count = region_chunks.size();
    int current_region = 0;
    // std::set<SVCall> combined_sv_calls;
    for (const auto& sub_region : region_chunks) {
        std::tuple<std::set<SVCall>, PrimaryMap, SuppMap> region_data = this->detectCIGARSVs(fp_in, idx, bamHdr, sub_region);
        std::set<SVCall>& subregion_sv_calls = std::get<0>(region_data);
        PrimaryMap& primary_map = std::get<1>(region_data);
        SuppMap& supp_map = std::get<2>(region_data);
        // std::cout << "Merge CIGAR SV calls from " << sub_region << "..." << std::endl;
        mergeSVs(subregion_sv_calls);
        int region_sv_count = getSVCount(subregion_sv_calls);
        // printMessage("Total SVs detected from CIGAR string: " + std::to_string(region_sv_count));

        // Run copy number variant predictions on the SVs detected from the
        // CIGAR string, using a minimum CNV length threshold
        if (region_sv_count > 0) {
            // std::cout << "Running copy number variant detection from CIGAR string SVs..." << std::endl;
            cnv_caller.runCIGARCopyNumberPrediction(chr, subregion_sv_calls, min_cnv_length, hmm);
        }

        // Run split-read SV and copy number variant predictions
        // std::cout << "Detecting copy number variants from split reads..." << std::endl;
        this->detectSVsFromSplitReads(subregion_sv_calls, primary_map, supp_map, cnv_caller, hmm);

        // Merge the SV calls from the current region
        // std::cout << "Merge SV calls from " << sub_region << "..." << std::endl;
        mergeSVs(subregion_sv_calls);

        // Combine the SV calls from the current region
        // std::cout << "Combining SV calls from " << sub_region << "..." << std::endl;
        concatenateSVCalls(combined_sv_calls, subregion_sv_calls);
        current_region++;
        printMessage("Completed " + std::to_string(current_region) + " of " + std::to_string(region_count) + " region(s) for chromosome " + chr + "...");
    }

    // Clean up the BAM file, header, and index
    hts_idx_destroy(idx);
    bam_hdr_destroy(bamHdr);
    sam_close(fp_in);
}

std::unordered_map<std::string, std::set<SVCall>> SVCaller::run()
{
    // Get the chromosomes to process
    std::vector<std::string> chromosomes;
    if (this->input_data.getChromosome() != "") {
        chromosomes.push_back(this->input_data.getChromosome());
    } else {
        chromosomes = this->input_data.getRefGenomeChromosomes();
    }

    // Ignore all alternate contigs (contains 'alt', 'GL', 'NC', 'hs', etc.)
    chromosomes.erase(std::remove_if(chromosomes.begin(), chromosomes.end(), [](const std::string& chr) {
        return chr.find("alt") != std::string::npos || chr.find("GL") != std::string::npos || chr.find("NC") != std::string::npos || chr.find("hs") != std::string::npos;
    }), chromosomes.end());

    // Read the HMM from the file
    std::string hmm_filepath = this->input_data.getHMMFilepath();
    std::cout << "Reading HMM from file: " << hmm_filepath << std::endl;
    CHMM hmm = ReadCHMM(hmm_filepath.c_str());

    // Set up threads for processing each chromosome
    std::vector<std::future<void>> futures;
    std::unordered_map<std::string, std::set<SVCall>> whole_genome_sv_calls;
    std::mutex sv_mutex;

    // Set a thread count for processing each chromosome. Keep it low to avoid
    // memory issues.
    int max_threads = 6;  // Number of chromosomes to process in parallel
    int batch_count = 0;
    int completed_threads = 0;
    int chr_count = chromosomes.size();
    for (const auto& chr : chromosomes) {
        printMessage("Launching thread for chromosome " + chr + "...");
        futures.push_back(std::async(std::launch::async, [&]() {
            std::set<SVCall> sv_calls;
            this->processChromosome(chr, this->input_data.getLongReadBam(), hmm, sv_calls, this->input_data.getMinCNVLength());
            {
                std::lock_guard<std::mutex> lock(sv_mutex);
                whole_genome_sv_calls[chr] = std::move(sv_calls);
            }
        }
        ));
        batch_count++;
        if (batch_count >= max_threads || batch_count >= chr_count) {
            // Wait for all threads to finish
            // printMessage("Waiting for all threads to finish for " + std::to_string(batch_count) + " chromosome(s)...");
            for (auto& future : futures) {
                future.get();
                completed_threads++;
                printMessage("Completed " + std::to_string(completed_threads) + " of " + std::to_string(chr_count) + " chromosome(s)");
            }
            // completed_threads += batch_count;
            // printMessage("Completed " + std::to_string(completed_threads) + " of " + std::to_string(chr_count) + " chromosome(s)");
            batch_count = 0;
            futures.clear();
        }
    }

    // Wait for remaining threads to finish
    if (futures.size() > 0) {
        // printMessage("Waiting for remaining threads to finish for " + std::to_string(futures.size()) + " chromosome(s)...");
        for (auto& future : futures) {
            future.get();
            completed_threads++;
            printMessage("[TEST] Completed " + std::to_string(completed_threads) + " of " + std::to_string(chr_count) + " chromosome(s)");
        }
        // completed_threads += futures.size();
        // printMessage("Completed " + std::to_string(completed_threads) + " of " + std::to_string(chr_count) + " chromosome(s)");
    }

    // // Loop through each region and detect SVs in chunks
    // std::string bam_filepath = this->input_data.getLongReadBam();
    // int chr_count = chromosomes.size();
    // std::cout << "Detecting SVs from " << chr_count << " chromosome(s)..." << std::endl;
    // // int thread_count = this->input_data.getThreadCount();
    // int thread_count = chr_count;
    // int min_cnv_length = this->input_data.getMinCNVLength();
    // for (const auto& chr : chromosomes) {
    //     printMessage("Launching thread for chromosome " + chr + "...");
    //     futures.push_back(std::async(std::launch::async, [&]() {
    //         std::set<SVCall> sv_calls;
    //         this->processChromosome(chr, bam_filepath, hmm, sv_calls, min_cnv_length);
    //         {
    //             std::lock_guard<std::mutex> lock(sv_mutex);
    //             whole_genome_sv_calls[chr] = std::move(sv_calls);
    //         }
    //     }
    //     ));
    // }

    // // Wait for all threads to finish
    // printMessage("Waiting for all threads to finish for " + std::to_string(chr_count) + " chromosome(s)...");
    // int threads_finished = 0;
    // for (auto& future : futures) {
    //     try{
    //         // future.wait();
    //         future.get();  // Wait and handle exceptions
    //         threads_finished++;
    //         printMessage("Completed " + std::to_string(threads_finished) + " of " + std::to_string(thread_count) + " threads...");
    //     } catch (const std::exception& e) {
    //         std::cerr << "Error in thread: " << e.what() << std::endl;
    //     }
    // }
    printMessage("All threads have finished.");

    // Print the total number of SVs detected for each chromosome
    uint32_t total_sv_count = 0;
    for (const auto& entry : whole_genome_sv_calls) {
        std::string chr = entry.first;
        int sv_count = getSVCount(entry.second);
        total_sv_count += sv_count;
        printMessage("Total SVs detected for chromosome " + chr + ": " + std::to_string(sv_count));
    }
    printMessage("Total SVs detected for all chromosomes: " + std::to_string(total_sv_count));

    // Save to VCF
    std::cout << "Saving SVs to VCF..." << std::endl;
    this->saveToVCF(whole_genome_sv_calls);

    return whole_genome_sv_calls;
}


// Detect SVs from split read alignments
void SVCaller::detectSVsFromSplitReads(std::set<SVCall>& sv_calls, PrimaryMap& primary_map, SuppMap& supp_map, CNVCaller& cnv_caller, CHMM hmm)
{
    // Find split-read SV evidence
    int sv_count = 0;
    int min_cnv_length = this->input_data.getMinCNVLength();
    for (const auto& entry : primary_map) {
        std::string qname = entry.first;
        AlignmentData primary_alignment = entry.second;
        std::string primary_chr = std::get<0>(primary_alignment);
        int32_t primary_start = std::get<1>(primary_alignment);
        int32_t primary_end = std::get<2>(primary_alignment);
        std::unordered_map<int, int> primary_match_map = std::get<6>(primary_alignment);

        // Skip primary alignments that do not have supplementary alignments
        if (supp_map.find(qname) == supp_map.end()) {
            continue;
        }

        // Find the largest supplementary alignment, and also identify inversions
        AlignmentData largest_supp_alignment = supp_map[qname][0];
        int32_t largest_supp_length = 0;
        for (auto it = supp_map[qname].begin(); it != supp_map[qname].end(); ++it) {
            const auto& supp_chr = std::get<0>(*it);
            if (primary_chr != supp_chr) {
                continue;  // Skip supplementary alignments on different chromosomes
            }
            int32_t supp_start = std::get<1>(*it);
            int32_t supp_end = std::get<2>(*it);
            int32_t supp_length = supp_end - supp_start + 1;
            if (supp_length > largest_supp_length) {
                largest_supp_length = supp_length;
                largest_supp_alignment = *it;
            }

            // Inversion detection
            bool is_opposite_strand = std::get<7>(primary_alignment) != std::get<7>(*it);
            if (is_opposite_strand) {
                if (supp_length >= min_cnv_length) {
                    SVCandidate sv_candidate(supp_start+1, supp_end+1, ".");
                    std::tuple<double, SVType, std::string, bool> result = cnv_caller.runCopyNumberPrediction(primary_chr, sv_candidate, hmm);
                    double supp_lh = std::get<0>(result);
                    SVType supp_type = std::get<1>(result);
                    if (supp_type == SVType::NEUTRAL) {
                        addSVCall(sv_calls, (uint32_t)(supp_start+1), (uint32_t)(supp_end+1), "INV", ".", "HMM", "./.", supp_lh);
                        sv_count++;
                    } else if (supp_type == SVType::DUP) {
                        addSVCall(sv_calls, (uint32_t)(supp_start+1), (uint32_t)(supp_end+1), "INVDUP", ".", "HMM", "./.", supp_lh);
                        sv_count++;
                    }
                } else {
                    // Add the inversion without running copy number predictions
                    // (too small for predictions)
                    addSVCall(sv_calls, (uint32_t)(supp_start+1), (uint32_t)(supp_end+1), "INV", ".", "REV", "./.", 0.0);
                    sv_count++;
                }
            }
        }

        // Trim overlapping alignments
        int32_t supp_start = std::get<1>(largest_supp_alignment);
        int32_t supp_end = std::get<2>(largest_supp_alignment);
        bool primary_before_supp = primary_start < supp_start;
        trimOverlappingAlignments(primary_alignment, largest_supp_alignment);

        // Create the SV candidate using both alignments
        supp_start = std::get<1>(largest_supp_alignment);
        supp_end = std::get<2>(largest_supp_alignment);
        primary_start = std::get<1>(primary_alignment);
        primary_end = std::get<2>(primary_alignment);
        SVCandidate split_boundary;
        SVCandidate split_gap;
        bool gap_exists = false;
        int32_t boundary_left, boundary_right, gap_left, gap_right;
        if (primary_before_supp) {
            boundary_left = primary_start+1;
            boundary_right = supp_end+1;
            gap_left = primary_end+1;
            gap_right = supp_start+1;
            gap_exists = primary_end < supp_start;
        } else {
            boundary_left = supp_start+1;
            boundary_right = primary_end+1;
            gap_left = supp_end+1;
            gap_right = primary_start+1;
            gap_exists = supp_end < primary_start;
        }
        
        // Run copy number variant predictions on the boundary if large enough
        if (boundary_right - boundary_left >= min_cnv_length) {
            split_boundary = SVCandidate(boundary_left, boundary_right, ".");
            std::tuple<double, SVType, std::string, bool> bd_result = cnv_caller.runCopyNumberPrediction(primary_chr, split_boundary, hmm);
            double bd_lh = std::get<0>(bd_result);
            SVType bd_type = std::get<1>(bd_result);

            // Run copy number variant predictions on the gap if it exists
            if (gap_exists && gap_right - gap_left >= min_cnv_length) {
                split_gap = SVCandidate(gap_left, gap_right, ".");
                std::tuple<double, SVType, std::string, bool> gap_result = cnv_caller.runCopyNumberPrediction(primary_chr, split_gap, hmm);
                double gap_lh = std::get<0>(gap_result);
                SVType gap_type = std::get<1>(gap_result);

                // If higher likelihood than the boundary, add the gap as the SV call
                if (gap_lh > bd_lh) {
                    addSVCall(sv_calls, (uint32_t)(gap_left), (uint32_t)(gap_right), getSVTypeString(gap_type), ".", "GAP", "./.", gap_lh);
                    sv_count++;
                } else {
                    // Add the boundary as the SV call
                    addSVCall(sv_calls, (uint32_t)(boundary_left), (uint32_t)(boundary_right), getSVTypeString(bd_type), ".", "BOUNDARY", "./.", bd_lh);
                    sv_count++;
                }
            } else {
                // Add the boundary as the SV call
                addSVCall(sv_calls, (uint32_t)(boundary_left), (uint32_t)(boundary_right), getSVTypeString(bd_type), ".", "BOUNDARY", "./.", bd_lh);
                sv_count++;
            }
        }
    }

    // Print the number of SVs detected from split-read alignments
    if (sv_count > 0) {
        std::cout << "Found " << sv_count << " SVs from split-read alignments" << std::endl;
    }
}

void SVCaller::saveToVCF(const std::unordered_map<std::string, std::set<SVCall> >& sv_calls)
{
    std::cout << "Creating VCF writer..." << std::endl;
    // std::string output_vcf = output_dir + "/output.vcf";
    std::string output_vcf = this->input_data.getOutputDir() + "/output.vcf";
    std::cout << "Writing VCF file to " << output_vcf << std::endl;
	std::ofstream vcf_stream(output_vcf);
    if (!vcf_stream.is_open()) {
        throw std::runtime_error("Failed to open VCF file for writing.");
    }
    std::string sample_name = "SAMPLE";

    std::cout << "Getting reference genome filepath..." << std::endl;
    try {
        std::string ref_fp = this->input_data.getRefGenome().getFilepath();
        std::cout << "Reference genome filepath: " << ref_fp << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return;
    }

    // Set the header lines
    std::cout << "Getting reference genome header..." << std::endl;
    const std::string contig_header = this->input_data.getRefGenome().getContigHeader();
    std::vector<std::string> header_lines = {
        std::string("##reference=") + 
        contig_header,
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">",
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">",
        "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">",
        "##INFO=<ID=SVMETHOD,Number=1,Type=String,Description=\"Method used to call the structural variant\">",
        "##INFO=<ID=ALN,Number=1,Type=String,Description=\"Alignment type used to call the structural variant\">",
        "##INFO=<ID=CLIPSUP,Number=1,Type=Integer,Description=\"Clipped base support at the start and end positions\">",
        "##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description=\"Number of reads supporting the variant\">",
        "##INFO=<ID=REPTYPE,Number=1,Type=String,Description=\"Repeat type\">",
        "##INFO=<ID=HMM,Number=1,Type=Float,Description=\"HMM likelihood\">",
        "##FILTER=<ID=PASS,Description=\"All filters passed\">",
        "##FILTER=<ID=LowQual,Description=\"Low quality\">",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">"
    };

    std::cout << "Writing VCF header..." << std::endl;

    // Add the file format
    std::string file_format = "##fileformat=VCFv4.2";
    vcf_stream << file_format << std::endl;

    // Add date and time
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];
    time (&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(buffer, sizeof(buffer), "%Y%m%d", timeinfo);
    vcf_stream << "##fileDate=" << buffer << std::endl;

    // Add source
    std::string source = "##source=ContexSV";
    vcf_stream << source << std::endl;

    // Loop over the header metadata lines
    for (const auto &line : header_lines) {
        vcf_stream << line << std::endl;
    }

    // Add the header line
    std::string header_line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE";
    vcf_stream << header_line << std::endl;

    // Flush the stream to ensure that the header is written
    //this->file_stream.flush();

    std::cout << "Saving SV calls to " << output_vcf << std::endl;
    std::string sv_method = "CONTEXTSVv0.1";
    int skip_count = 0;
    int total_count = 0;
    for (const auto& pair : sv_calls) {
        std::string chr = pair.first;
        const std::set<SVCall>& sv_calls = pair.second;
        std::cout << "Saving SV calls for " << chr << "..." << std::endl;
        for (const auto& sv_call : sv_calls) {
            // Get the SV candidate and SV info
            uint32_t start = sv_call.start;
            uint32_t end = sv_call.end;
            std::string sv_type_str = sv_call.sv_type;
            std::string genotype = sv_call.genotype;
            std::string data_type_str = sv_call.data_type;
            std::string alt_allele = sv_call.alt_allele;
            double hmm_likelihood = sv_call.hmm_likelihood;
            int sv_length = end - start;
            if (sv_type_str == "DEL") {
            	sv_length++;
        	}
            int read_support = sv_call.support;
            int read_depth = 0;
            // SVType sv_type = sv_call.sv_type;
            // SVCandidate candidate = sv_call.first;
            // SVInfo info = sv_call.second;
            // SVType sv_type = info.sv_type;
            // int read_support = info.read_support;
            // int read_depth = info.read_depth;
            // int read_depth = 0;
            // int read_support = 0;
            // int sv_length = info.sv_length;
            // std::set<std::string> data_type = info.data_type;
            // std::string genotype = info.genotype;
            // double hmm_likelihood = info.hmm_likelihood;

            // Convert the data type set to a string
            // std::string data_type_str = "";
            // for (auto const& type : data_type) {
            //     data_type_str += type + ",";
            // }

            // Get the CHROM, POS, END, and ALT
            // uint32_t pos = std::get<0>(candidate);
            // uint32_t end = std::get<1>(candidate);

            // If the SV type is unknown, skip it
            if (sv_type_str == "UNKNOWN" || sv_type_str == "NEUTRAL") {
                skip_count += 1;
                continue;
            } else {
                total_count += 1;
            }

            // Process by SV type
            std::string ref_allele = ".";
            // std::string alt_allele = ".";
            std::string repeat_type = "NA";

            // Deletion
            if (sv_type_str == "DEL") {
                // Get the deleted sequence from the reference genome, also including the preceding base
                int64_t preceding_pos = (int64_t) std::max(1, (int) start-1);  // Make sure the position is not negative
                // ref_allele = ref_genome.query(chr, preceding_pos, end);
                ref_allele = this->input_data.queryRefGenome(chr, preceding_pos, end);

                // Use the preceding base as the alternate allele 
                if (ref_allele != "") {
                    alt_allele = ref_allele.at(0);
                } else {
                    alt_allele = "<DEL>";  // Symbolic allele
                    std::cerr << "Warning: Reference allele is empty for deletion at " << chr << ":" << start << "-" << end << std::endl;
                }

                sv_length = -1 * sv_length;  // Negative length for deletions

                start = preceding_pos;  // Update the position to the preceding base

            // Other types (duplications, insertions, inversions)
            } else {
                // Use the preceding base as the reference allele
                int64_t preceding_pos = (int64_t) std::max(1, (int) start-1);  // Make sure the position is not negative
                // ref_allele = ref_genome.query(chr, preceding_pos,
                // preceding_pos);
                ref_allele = this->input_data.queryRefGenome(chr, preceding_pos, preceding_pos);

                // Format novel insertions
                if (sv_type_str == "INS") {
                    // Check if in symbolic form
                    if (alt_allele != "<INS>") {
                        // Use the insertion sequence as the alternate allele
                        // alt_allele = std::get<2>(candidate);
                        alt_allele.insert(0, ref_allele);
                    }
                    start = preceding_pos;  // Update the position to the preceding base

                    // Update the end position to the start position to change from
                    // query to reference coordinates for insertions
                    end = start;
                }
            }

            // Create the VCF parameter strings
            // int clipped_base_support = this->getClippedBaseSupport(chr, pos,
            // end);
            int clipped_base_support = 0;
            // std::string sv_type_str = getSVTypeString(sv_type);
            std::string info_str = "END=" + std::to_string(end) + ";SVTYPE=" + sv_type_str + \
                ";SVLEN=" + std::to_string(sv_length) + ";SUPPORT=" + std::to_string(read_support) + \
                ";SVMETHOD=" + sv_method + ";ALN=" + data_type_str + ";CLIPSUP=" + std::to_string(clipped_base_support) + \
                ";REPTYPE=" + repeat_type + ";HMM=" + std::to_string(hmm_likelihood);
                
            std::string format_str = "GT:DP";
            std::string sample_str = genotype + ":" + std::to_string(read_depth);
            std::vector<std::string> samples = {sample_str};

            // Write the SV call to the file (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLES)
            vcf_stream << chr << "\t" << start << "\t" << "." << "\t" << ref_allele << "\t" << alt_allele << "\t" << "." << "\t" << "PASS" << "\t" << info_str << "\t" << format_str << "\t" << samples[0] << std::endl;
            if (total_count % 1000 == 0)
            {
            	std::cout << "Wrote SV at " << chr << ": " << start << ", total=" << total_count << std::endl;
        	}
        }
    }

    // Print the number of SV calls skipped
    std::cout << "Finished writing VCF file. Total SV calls: " << total_count << ", skipped: " << skip_count << " with unknown SV type" << std::endl;
}

void SVCaller::trimOverlappingAlignments(AlignmentData& primary_alignment, AlignmentData& supp_alignment)
{
    // Get the start and end read positions for the primary and supplementary
    // alignments
    int32_t primary_query_start = std::get<4>(primary_alignment);
    int32_t primary_query_end = std::get<5>(primary_alignment);
    int32_t supp_query_start = std::get<4>(supp_alignment);
    int32_t supp_query_end = std::get<5>(supp_alignment);
    std::unordered_map<int, int>& primary_match_map = std::get<6>(primary_alignment);
    std::unordered_map<int, int>& supp_match_map = std::get<6>(supp_alignment);
    int32_t primary_alignment_start = std::get<1>(primary_alignment);
    int32_t primary_alignment_end = std::get<2>(primary_alignment);
    int32_t supp_alignment_start = std::get<1>(supp_alignment);
    int32_t supp_alignment_end = std::get<2>(supp_alignment);

    // Check if the alignments overlap
    bool primary_before_supp = primary_query_start < supp_query_start;
    if (primary_before_supp) {
        // Primary before supplementary in the query
        if (primary_query_end >= supp_query_start) {
            // Calculate the mismatch rates at the overlapping region
            double primary_mismatch_rate = this->calculateMismatchRate(primary_match_map, supp_query_start, primary_query_end);
            double supp_mismatch_rate = this->calculateMismatchRate(supp_match_map, supp_query_start, primary_query_end);
            int32_t overlap_length = primary_query_end - supp_query_start + 1;

            // Trim the ailgnment with the higher mismatch rate
            if (primary_mismatch_rate > supp_mismatch_rate) {
                // Trim the end of the primary alignment
                std::get<2>(primary_alignment) = primary_alignment_end - overlap_length;
            } else {
                // Trim the beginning of the supplementary alignment
                std::get<1>(supp_alignment) = supp_alignment_start + overlap_length;
            }
        }
    } else {
        // Supplementary before primary in the query
        if (supp_query_end >= primary_query_start) {
            // Calculate the mismatch rates at the overlapping region
            double primary_mismatch_rate = this->calculateMismatchRate(primary_match_map, primary_query_start, supp_query_end);
            double supp_mismatch_rate = this->calculateMismatchRate(supp_match_map, primary_query_start, supp_query_end);
            int32_t overlap_length = supp_query_end - primary_query_start + 1;

            // Trim the ailgnment with the higher mismatch rate
            if (supp_mismatch_rate > primary_mismatch_rate) {
                // Trim the end of the supplementary alignment
                std::get<2>(supp_alignment) = supp_alignment_end - overlap_length;
            } else {
                // Trim the beginning of the primary alignment
                std::get<1>(primary_alignment) = primary_alignment_start + overlap_length;
            }
        }
    }
}
