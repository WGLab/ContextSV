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
    // // Open the BAM file
    // std::string bam_filepath = this->input_data.getLongReadBam();
    // samFile *fp_in = sam_open(bam_filepath.c_str(), "r");
    // if (fp_in == NULL) {
    //     std::cerr << "ERROR: failed to open " << bam_filepath << std::endl;
    //     exit(1);
    // }

    // // Load the header for the BAM file
    // bam_hdr_t *bamHdr = sam_hdr_read(fp_in);
    // if (!bamHdr) {
    //     sam_close(fp_in);
    //     throw std::runtime_error("ERROR: failed to read header for " + bam_filepath);
    // }

    // // Load the index for the BAM file
    // hts_idx_t *idx = sam_index_load(fp_in, bam_filepath.c_str());
    // if (!idx) {
    //     bam_hdr_destroy(bamHdr);
    //     sam_close(fp_in);
    //     throw std::runtime_error("ERROR: failed to load index for " + bam_filepath);
    // }

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
    // SVData sv_calls;
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
    
    // hts_itr_destroy(itr);
    // bam_destroy1(bam1);
    // hts_idx_destroy(idx);
    // bam_hdr_destroy(bamHdr);
    // sam_close(fp_in);

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

                // Loop from the leftmost position of the insertion (pos-op_len)
                // to the rightmost position of the insertion (pos+op_len-1) and
                // calculate the sequence identity at each window of the
                // insertion length to identify potential duplications.

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
                    // sv_calls.add(chr, ref_pos, ref_end, SVType::DUP,
                    // ins_seq_str, "CIGARDUP", "./.", default_lh);
                    //printMessage("[TEST] FOUND CIGAR DUP");
                    // sv_calls.insert(SVCall{(uint32_t)ref_pos,
                    // (uint32_t)ref_end, "DUP", ins_seq_str, "CIGARDUP", "./.",
                    // default_lh});
                    addSVCall(sv_calls, (uint32_t)ref_pos, (uint32_t)ref_end, "DUP", ins_seq_str, "CIGARDUP", "./.", default_lh);
                } else {
                    // sv_calls.add(chr, ref_pos, ref_end, SVType::INS, ins_seq_str, "CIGARINS", "./.", default_lh);
                    // sv_calls.insert(SVCall{(uint32_t)ref_pos,
                    // (uint32_t)ref_end, "INS", ins_seq_str, "CIGARINS", "./.",
                    // default_lh});
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
                // sv_calls.add(chr, ref_pos, ref_end, SVType::DEL, ".",
                // "CIGARDEL", "./.", default_lh);  // Add to SV calls (1-based)
                // sv_calls.insert(SVCall{(uint32_t)ref_pos, (uint32_t)ref_end,
                // "DEL", ".", "CIGARDEL", "./.", default_lh});
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
                std::cerr << "ERROR: Sequence lengths do not match" << std::endl;
                exit(1);
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
            std::cerr << "ERROR: Unknown CIGAR operation " << op << std::endl;
            exit(1);
        }

        // Update the query position based on the CIGAR operation (M, I, S, H)
        if (op == BAM_CMATCH || op == BAM_CINS || op == BAM_CSOFT_CLIP || op == BAM_CEQUAL || op == BAM_CDIFF) {
            query_pos += op_len;
        } else if (op == BAM_CDEL || op == BAM_CREF_SKIP || op == BAM_CHARD_CLIP || op == BAM_CPAD) {
            // Do nothing
        } else {
            std::cerr << "ERROR: Unknown CIGAR operation " << op << std::endl;
            exit(1);
        }
    }

    query_end = query_pos;  // Last alignment position in the query

    return std::tuple<std::unordered_map<int, int>, int32_t, int32_t>(query_match_map, query_start, query_end);
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

    // [TEST] Only process the last N chromosomes
    // int last_n = 3;
    // chromosomes = std::vector<std::string>(chromosomes.end()-last_n, chromosomes.end());
    // std::cout << "[DEBUG] Running last " << last_n << " chromosomes" << std::endl;
    // //chromosomes = std::vector<std::string>(chromosomes.end()-3, chromosomes.end());

    // Open the BAM file
    std::string bam_filepath = this->input_data.getLongReadBam();
    samFile *fp_in = sam_open(bam_filepath.c_str(), "r");
    if (!fp_in) {
        throw std::runtime_error("ERROR: failed to open " + bam_filepath);
    }

    // Load the header for the BAM file
    bam_hdr_t *bamHdr = sam_hdr_read(fp_in);
    if (!bamHdr) {
        sam_close(fp_in);
        throw std::runtime_error("ERROR: failed to read header for " + bam_filepath);
    }

    // Load the index for the BAM file
    hts_idx_t *idx = sam_index_load(fp_in, bam_filepath.c_str());
    if (!idx) {
        bam_hdr_destroy(bamHdr);
        sam_close(fp_in);
        throw std::runtime_error("ERROR: failed to load index for " + bam_filepath);
    }

    // Loop through each region and detect SVs in chunks
    int chr_count = chromosomes.size();
    int current_chr = 0;
    std::cout << "Detecting SVs from " << chr_count << " chromosome(s)..." << std::endl;
    int chunk_count = 100;  // Number of chunks to split the chromosome into
    // SVData sv_calls;
    // std::vector<std::map<SVCandidate, SVInfo>> sv_calls;
    // std::unordered_map<std::string, std::map<uint32_t, uint32_t>> sv_calls;
    uint32_t total_sv_count = 0;
    std::unordered_map<std::string, std::set<SVCall>> whole_genome_sv_calls;
    int min_cnv_length = this->input_data.getMinCNVLength();
    for (const auto& chr : chromosomes) {
        std::cout << "Running SV detection for chromosome " << chr << "..." << std::endl;

        // Split the chromosome into chunks
        std::vector<std::string> region_chunks;
        if (this->input_data.isRegionSet()) {

            // Use one chunk for the specified region
            std::pair<int32_t, int32_t> region = this->input_data.getRegion();
            int region_start = region.first;
            int region_end = region.second;
            std::string chunk = chr + ":" + std::to_string(region_start) + "-" + std::to_string(region_end);
            region_chunks.push_back(chunk);
            std::cout << "Using specified region " << chunk << "..." << std::endl;
            
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
            std::cout << "Split chromosome " << chr << " into " << region_chunks.size() << " chunks of size " << chunk_size << "..." << std::endl;
        }

        // Load chromosome data for copy number predictions
        std::cout << "Loading chromosome data for copy number predictions..." << std::endl;
        CNVCaller cnv_caller(this->input_data);
        cnv_caller.loadChromosomeData(chr);

        // Process each chunk one at a time
        std::cout << "Processing " << region_chunks.size() << " region(s) for chromosome " << chr << "..." << std::endl;
        int region_count = region_chunks.size();
        int current_region = 0;
        std::set<SVCall> combined_sv_calls;
        for (const auto& sub_region : region_chunks) {
            // std::cout << "Detecting CIGAR string SVs from " << sub_region << "..." << std::endl;
            // std::tuple<std::set<SVCall>, PrimaryMap, SuppMap> region_data = this->detectCIGARSVs(sub_region);
            std::tuple<std::set<SVCall>, PrimaryMap, SuppMap> region_data = this->detectCIGARSVs(fp_in, idx, bamHdr, sub_region);
            std::set<SVCall>& subregion_sv_calls = std::get<0>(region_data);
            PrimaryMap& primary_map = std::get<1>(region_data);
            SuppMap& supp_map = std::get<2>(region_data);
            // SVData& subregion_sv_calls = std::get<0>(region_data);
            // PrimaryMap& primary_map = std::get<1>(region_data);
            // SuppMap& supp_map = std::get<2>(region_data);
            // int region_sv_count = subregion_sv_calls.totalCalls();
            // if (region_sv_count > 0) {
            //     std::cout << "Detected " << region_sv_count << " CIGAR SVs from " << sub_region << "..." << std::endl;
            // }
            // int region_sv_count = subregion_sv_calls.count();
            int region_sv_count = getSVCount(subregion_sv_calls);
            printMessage("Total SVs detected from CIGAR string: " + std::to_string(region_sv_count));

            // Run copy number variant predictions on the SVs detected from the
            // CIGAR string, using a minimum CNV length threshold
            // std::cout << "Detecting copy number variants from CIGAR string SVs..." << std::endl;
            // std::map<SVCandidate, SVInfo>& cigar_svs = subregion_sv_calls.getChromosomeSVs(chr);
            // if (cigar_svs.size() > 0) {
            if (region_sv_count > 0) {
                std::cout << "Running copy number variant detection from CIGAR string SVs..." << std::endl;
                // cnv_caller.runCIGARCopyNumberPrediction(chr, cigar_svs,
                // min_cnv_length);
                cnv_caller.runCIGARCopyNumberPrediction(chr, subregion_sv_calls, min_cnv_length);
            }

            // Run split-read SV and copy number variant predictions
            std::cout << "Detecting copy number variants from split reads..." << std::endl;
            this->detectSVsFromSplitReads(subregion_sv_calls, primary_map, supp_map, cnv_caller);
            // sv_calls.concatenate(subregion_sv_calls);  // Add the calls to the
            // main set
            // sv_calls.emplace_back(subregion_sv_calls);

            // Merge the SV calls from the current region
            std::cout << "Merge SV calls from " << sub_region << "..." << std::endl;
            mergeSVs(subregion_sv_calls);

            // Combine the SV calls from the current region
            std::cout << "Combining SV calls from " << sub_region << "..." << std::endl;
            concatenateSVCalls(combined_sv_calls, subregion_sv_calls);
            std::cout << "Completed " << ++current_region << " of " << region_count << " region(s)..." << std::endl;

            // [TEST] Break after the first region
            // std::cout << "[DEBUG] Breaking after the first region" << std::endl;
            // break;
        }

        std::cout << "Completed " << ++current_chr << " of " << chr_count << " chromosome(s)..." << std::endl;
        int chr_sv_count = getSVCount(combined_sv_calls);
        whole_genome_sv_calls[chr] = combined_sv_calls;
        std::cout << "Total SVs detected for chromosome " << chr << ": " << chr_sv_count << std::endl;
        total_sv_count += chr_sv_count;
        std::cout << "Cumulative total SVs: " << total_sv_count << std::endl;
        // std::cout << "Completed " << region_count << " of " << chr_count << " chromosome(s)" << std::endl;
    }

    // Clean up the BAM file, header, and index
    hts_idx_destroy(idx);
    bam_hdr_destroy(bamHdr);
    sam_close(fp_in);

    // SVData sv_calls_combined;
    // for (const auto& subregion_sv_calls : sv_calls) {
    //     sv_calls_combined.concatenate(subregion_sv_calls);
    // }

    // Save to VCF
    std::cout << "Saving SVs to VCF..." << std::endl;
    this->saveToVCF(whole_genome_sv_calls);

    return whole_genome_sv_calls;
}


// Detect SVs from split read alignments
void SVCaller::detectSVsFromSplitReads(std::set<SVCall>& sv_calls, PrimaryMap& primary_map, SuppMap& supp_map, CNVCaller& cnv_caller)
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
        int32_t primary_query_start = std::get<4>(primary_alignment);
        int32_t primary_query_end = std::get<5>(primary_alignment);
        std::unordered_map<int, int> primary_match_map = std::get<6>(primary_alignment);
        // bool primary_strand = std::get<7>(primary_alignment);
        if (supp_map.find(qname) == supp_map.end()) {
            continue;
        }

        // // Resolve overlaps between the primary and supplementary query
        // // sequences
        // for (auto it = supp_map[qname].begin(); it != supp_map[qname].end(); ++it) {
        //     std::string supp_chr = std::get<0>(*it);
        //     // int32_t supp_start = std::get<1>(*it);
        //     // int32_t supp_end = std::get<2>(*it);
        //     int32_t supp_query_start = std::get<4>(*it);
        //     int32_t supp_query_end = std::get<5>(*it);
        //     std::unordered_map<int, int> supp_match_map = std::get<6>(*it);
        //     // bool supp_strand = std::get<7>(*it);

        //     // Resolve overlaps between the primary and supplementary query
        //     // sequences
        //     if (primary_query_start < supp_query_end && primary_query_end > supp_query_start || supp_query_start < primary_query_end && supp_query_end > primary_query_start) {

        //         // Calculate the mismatch rate for each alignment at the overlap
        //         double primary_mismatch_rate = this->calculateMismatchRate(primary_match_map, overlap_start, overlap_end-1);
        //         double supp_mismatch_rate = this->calculateMismatchRate(supp_match_map, overlap_start, overlap_end-1);
        //         // std::cout << "Primary mismatch rate: " << primary_mismatch_rate << std::endl;
        //         // std::cout << "Supplementary mismatch rate: " << supp_mismatch_rate << std::endl;

        //         // Trim the overlap from the alignment with the higher mismatch
        //         // rate
        //         if (primary_mismatch_rate > supp_mismatch_rate) {
        //             if (overlap_start == primary_query_start) {
        //                 primary_start += overlap_length;
        //             } else if (overlap_end == primary_query_end) {
        //                 primary_end -= overlap_length;
        //             }

        //         } else {
        //             if (overlap_start == supp_query_start) {
        //                 // supp_start += overlap_length;
        //                 // Update the value in the supp map
        //                 std::get<1>(*it) += overlap_length;
        //             } else if (overlap_end == supp_query_end) {
        //                 // supp_end -= overlap_length;
        //                 // Update the value in the supp map
        //                 std::get<2>(*it) -= overlap_length;
        //             }
        //         }
        //     }
        // }

        // Remove supplementary alignments that are not on the same chromosome
        // as the primary alignment
        for (auto it = supp_map[qname].begin(); it != supp_map[qname].end();) {
            if (std::get<0>(*it) != primary_chr) {
                it = supp_map[qname].erase(it);
            } else {
                ++it;
            }
        }

        // Run copy number variant predictions on the primary alignment
        SVType primary_type = SVType::UNKNOWN;
        double primary_lh = std::numeric_limits<double>::lowest();
        int32_t primary_lh_t = 0;
        if (primary_end - primary_start >= min_cnv_length) {
            SVCandidate sv_candidate(primary_start+1, primary_end+1, ".");
            // std::cout << "TEST5" << std::endl;
            std::tuple<double, SVType, std::string, bool> result = cnv_caller.runCopyNumberPrediction(primary_chr, sv_candidate);
            primary_lh = std::get<0>(result);
            // primary_log_likelihood /= (double)(primary_end - primary_start);  // Normalize the log likelihood by the length
            primary_type = std::get<1>(result);
        }

        // Loop through the supplementary alignments, find the largest
        // supplementary alignment, and the closest non-overlapping
        // supplementary alignment to the primary alignment
        AlignmentData largest_supp_alignment = supp_map[qname][0];
        AlignmentData closest_supp_alignment = supp_map[qname][0];
        int32_t largest_supp_length = 0;
        int32_t closest_supp_distance = std::numeric_limits<int32_t>::max();
        int32_t closest_supp_length = 0;
        for (auto it = supp_map[qname].begin(); it != supp_map[qname].end(); ++it) {
            // const auto& supp_chr = std::get<0>(*it);
            int32_t supp_start = std::get<1>(*it);
            int32_t supp_end = std::get<2>(*it);
            int32_t supp_length = supp_end - supp_start + 1;
            int32_t supp_distance = std::numeric_limits<int32_t>::max();
            if (supp_start > primary_end) {
                supp_distance = supp_start - primary_end;
            } else if (supp_end < primary_start) {
                supp_distance = primary_start - supp_end;
            }
            if (supp_length > largest_supp_length) {
                largest_supp_length = supp_length;
                largest_supp_alignment = *it;
            }
            if (supp_distance < closest_supp_distance) {
                closest_supp_length = supp_length;
                closest_supp_alignment = *it;
                closest_supp_distance = supp_distance;
            }
        }

        // Run copy number variant predictions on the largest supplementary
        // alignment
        double largest_supp_lh = std::numeric_limits<double>::lowest();
        SVType largest_supp_type = SVType::UNKNOWN;
        int largest_supp_lh_t = 0;
        if (largest_supp_length >= min_cnv_length) {
            SVCandidate sv_candidate(std::get<1>(largest_supp_alignment)+1, std::get<2>(largest_supp_alignment)+1, ".");
            // std::cout << "TEST1" << std::endl;
            std::tuple<double, SVType, std::string, bool> result = cnv_caller.runCopyNumberPrediction(primary_chr, sv_candidate);
            largest_supp_lh = std::get<0>(result);
            // largest_supp_log_likelihood /= (double)largest_supp_length;  // Normalize the log likelihood by the length
            largest_supp_type = std::get<1>(result);
        }

        // Run copy number variant predictions on the closest non-overlapping
        // supplementary alignment (if not the same as the largest)
        double closest_supp_lh = std::numeric_limits<double>::lowest();
        SVType closest_supp_type = SVType::UNKNOWN;
        int closest_supp_lh_t = 0;
        if (largest_supp_alignment != closest_supp_alignment) {
            if (closest_supp_length >= min_cnv_length) {
                SVCandidate sv_candidate(std::get<1>(closest_supp_alignment)+1, std::get<2>(closest_supp_alignment)+1, ".");
                // std::cout << "TEST2" << std::endl;
                std::tuple<double, SVType, std::string, bool> result = cnv_caller.runCopyNumberPrediction(primary_chr, sv_candidate);
                closest_supp_lh = std::get<0>(result);
                // closest_supp_log_likelihood /= (double)closest_supp_length;  // Normalize the log likelihood by the length
                closest_supp_type = std::get<1>(result);
                int32_t closest_supp_start = std::get<1>(closest_supp_alignment);
                int32_t closest_supp_end = std::get<2>(closest_supp_alignment);
            }
        }

        // Define constants representing read scenarios used for SV detection
        const int NOCALL = -1;  // Default
        const int PRIM_SUPP_BD = 0;  // Primary and supplementary boundary
        const int PRIM_SUPP_GAP = 1;  // Primary and supplementary gap
        const int SUPP_PRIM_BD = 2;  // Supplementary and primary boundary
        const int SUPP_PRIM_GAP = 3;  // Supplementary and primary gap

        // Loop through all the supplementary alignments and find the highest
        // likelihood prediction
        double best_split_aln_lh = std::numeric_limits<double>::lowest();
        double best_split_aln_lh_norm = std::numeric_limits<double>::lowest();
        // int best_split_aln_length = 0;
        SVType best_supp_type = SVType::UNKNOWN;
        std::pair<int32_t, int32_t> best_supp_candidate;
        AlignmentData& best_split_alignment = supp_map[qname][0];
        int best_scenario = NOCALL;
        for (auto it = supp_map[qname].begin(); it != supp_map[qname].end(); ++it) {
            int32_t supp_start = std::get<1>(*it);
            int32_t supp_end = std::get<2>(*it);
            bool primary_before_supp = primary_start < supp_start;

            // Create the SV candidate as the boundary of the primary and
            // supplementary alignments
            SVCandidate split_boundary;
            SVCandidate split_gap;
            bool invalid_gap = false;
            if (primary_before_supp) {
                split_boundary = SVCandidate(primary_start+1, supp_end+1, ".");

                // Check for an invalid gap (overlap)
                if (primary_end >= supp_start) {
                    invalid_gap = true;
                } else {
                    split_gap = SVCandidate(primary_end+1, supp_start+1, ".");
                }
                // split_gap = SVCandidate(primary_end+1, supp_start+1, ".");

            } else {
                split_boundary = SVCandidate(supp_start+1, primary_end+1, ".");

                // Check for an invalid gap (overlap)
                if (supp_end >= primary_start) {
                    invalid_gap = true;
                } else {
                    split_gap = SVCandidate(supp_end+1, primary_start+1, ".");
                }
            }

            // Create a vector of the two SV candidates, don't add the gap if
            // it is an overlap, or if either SV is less than the minimum CNV
            // length
            std::vector<SVCandidate> sv_candidates;
            if (!invalid_gap && std::get<1>(split_gap) - std::get<0>(split_gap) >= min_cnv_length) {
                sv_candidates.push_back(split_gap);
            }
            if (std::get<1>(split_boundary) - std::get<0>(split_boundary) >= min_cnv_length) {
                sv_candidates.push_back(split_boundary);
            }

            // Continue if no SV candidates
            if (sv_candidates.size() == 0) {
                continue;
            }

            // Run copy number variant predictions on both, and keep the
            // prediction with the highest normalized log likelihood
            double chosen_lh_norm = std::numeric_limits<double>::lowest();
            SVType chosen_type = SVType::UNKNOWN;
            std::pair<int32_t, int32_t> chosen_candidate;
            std::string chosen_candidate_str = "BOUNDARY";
            int split_scenario = NOCALL;
            for (const auto& sv_candidate : sv_candidates) {
            	// std::cout << "TEST3: primary = " << primary_start << ", " << primary_end << " supp = " << supp_start << ", " << supp_end << std::endl;
            	// std::cout << "Position: " << std::get<0>(sv_candidate) << ", " << std::get<1>(sv_candidate) << std::endl;
                std::tuple<double, SVType, std::string, bool> result = cnv_caller.runCopyNumberPrediction(primary_chr, sv_candidate);
                double current_lh = std::get<0>(result);
                SVType current_type = std::get<1>(result);

                // Normalize the log likelihood by the state sequence length
                double current_lh_norm = current_lh;// / (double)T;
                // if (sv_candidate == split_boundary) {
                //     std::cout << "Boundary candidate: " << current_lh_norm << ", before normalization: " << current_lh << std::endl;
                // } else if (sv_candidate == split_gap) {
                //     std::cout << "Gap candidate: " << current_lh_norm << ", before normalization: " << current_lh << std::endl;
                // }

                // Update the current SV candidate if the likelihood is higher
                if (current_type != SVType::UNKNOWN && current_lh_norm > chosen_lh_norm) {
                    chosen_lh_norm = current_lh_norm;
                    chosen_type = current_type;
                    chosen_candidate = std::make_pair(std::get<0>(sv_candidate), std::get<1>(sv_candidate));

                    // Update the candidate string
                    if (sv_candidate == split_boundary) {
                        chosen_candidate_str = "BOUNDARY";
                        if (primary_before_supp) {
                            split_scenario = PRIM_SUPP_BD;
                        } else {
                            split_scenario = SUPP_PRIM_BD;
                        }
                    } else if (sv_candidate == split_gap) {
                        chosen_candidate_str = "GAP";
                        if (primary_before_supp) {
                            split_scenario = PRIM_SUPP_GAP;
                        } else {
                            split_scenario = SUPP_PRIM_GAP;
                        }
                    }
                    // std::cout << "Updated candidate: " << chosen_candidate_str << " with likelihood: " << current_lh_norm << std::endl;
                } else if (current_type == SVType::UNKNOWN) {
                    // std::cerr << "ERROR: Unknown SV type" << std::endl;
                    // exit(1);
                }
            }

            // std::cout << "Chosen candidate: " << chosen_candidate_str << std::endl;

            // Continue if unknown SV type
            if (chosen_type == SVType::UNKNOWN) {
                // std::cerr << "ERROR: Unknown SV type" << std::endl;
                continue;
            }

            // If opposite strand, set the type to INV or INV_DUP
            bool same_strand = std::get<7>(*it) == std::get<7>(primary_alignment);
            if (!same_strand) {
                if (chosen_type == SVType::NEUTRAL) {
                    chosen_type = SVType::INV;
                } else if (chosen_type == SVType::DUP) {
                    chosen_type = SVType::INV_DUP;
                }
            }

            if (chosen_lh_norm > best_split_aln_lh_norm) {
                // best_supp_log_likelihood = supp_likelihood;
                // best_supp_log_likelihood /= (double)(sv_end - sv_start);  //
                // Normalize the log likelihood by the length
                // best_split_aln_lh = split_aln_lh;
                best_split_aln_lh_norm = chosen_lh_norm;
                // best_split_aln_length = split_aln_length;
                best_supp_type = chosen_type;
                best_supp_candidate = chosen_candidate;
                best_split_alignment = *it;
                best_scenario = split_scenario;
            } else if (chosen_lh_norm <= best_split_aln_lh_norm) {
                // std::cerr << "ERROR: split_aln_lh_norm is less than or equal to best_split_aln_lh_norm" << std::endl;
                // exit(1);
            }
        }

        // If the likelihood is equal to the lowest value, print an error
        if (best_split_aln_lh_norm == std::numeric_limits<double>::lowest()) {
            // std::cerr << "ERROR: best_supp_log_likelihood is the lowest value" << std::endl;
            // exit(1);
        }

        // Print the likelihoods
        // std::cout << "Primary log likelihood: " << primary_lh << std::endl;
        // std::cout << "Largest supplementary log likelihood: " << largest_supp_lh << std::endl;
        // std::cout << "Closest supplementary log likelihood: " << closest_supp_lh << std::endl;
        // // std::cout << "Best split alignment log likelihood: " << best_split_aln_lh << std::endl;
        // std::cout << "Best split alignment log likelihood (normalized): " << best_split_aln_lh_norm << std::endl;
        // std::cout << "Best scenario: " << best_scenario << std::endl;

        // Add the SV call with the highest likelihood prediction
        // 
        // Determine the normalized log likelihood for the combined alignments
        // by summing and normalizing the log likelihoods by the length
        double complex_lh = 0.0;
        double complex_lh_norm = 0.0;
        if (largest_supp_alignment == closest_supp_alignment) {
            int32_t complex_t = primary_lh_t + largest_supp_lh_t;
            complex_lh = primary_lh + largest_supp_lh;
            complex_lh_norm = complex_lh;// / complex_t;
        } else {
            int32_t complex_t = primary_lh_t + largest_supp_lh_t + closest_supp_lh_t;
            complex_lh = primary_lh + largest_supp_lh + closest_supp_lh;
            complex_lh_norm = complex_lh;// / complex_t;
        }
        // std::cout << "Complex log likelihood (normalized): " << complex_lh_norm << std::endl;

        // Compare the best split alignment likelihood to the complex likelihood
        // if (best_supp_log_likelihood > primary_log_likelihood || best_supp_log_likelihood > largest_supp_log_likelihood || best_supp_log_likelihood > closest_supp_log_likelihood) {
        if (best_split_aln_lh_norm > complex_lh_norm) {
            int32_t sv_start = best_supp_candidate.first;
            int32_t sv_end = best_supp_candidate.second;

            // Print an error and continue if the end is less than the start
            if (sv_end < sv_start) {
                std::cerr << "ERROR: SV end is less than the start: " << sv_start << " - " << sv_end << ", SV type: " << getSVTypeString(best_supp_type) << std::endl;
                continue;
            }

            // Resolve overlaps between the primary and supplementary query
            // sequences for deletions (not usually an issue for other types)
            if (best_supp_type == SVType::DEL) {
                AlignmentData& best_supp_alignment = best_split_alignment;
                int32_t supp_start = std::get<1>(best_supp_alignment);
                int32_t supp_end = std::get<2>(best_supp_alignment);
                int32_t supp_query_start = std::get<4>(best_supp_alignment);
                int32_t supp_query_end = std::get<5>(best_supp_alignment);
                std::unordered_map<int, int> supp_match_map = std::get<6>(best_supp_alignment);

                // Resolve overlaps between the primary and supplementary query
                // sequences
                // int32_t overlap_start = std::max(primary_query_start, supp_query_start);
                // int32_t overlap_end = std::min(primary_query_end, supp_query_end);
                // int32_t overlap_length = overlap_end - overlap_start;
                bool gap_present = primary_query_end < supp_query_start || supp_query_end < primary_query_start;
                if (!gap_present) {
                    int32_t overlap_start = std::max(primary_query_start, supp_query_start);
                    int32_t overlap_end = std::min(primary_query_end, supp_query_end);
                    int32_t overlap_length = overlap_end - overlap_start;

                    // Calculate the mismatch rate for each alignment at the overlap
                    double primary_mismatch_rate = this->calculateMismatchRate(primary_match_map, overlap_start, overlap_end);
                    double supp_mismatch_rate = this->calculateMismatchRate(supp_match_map, overlap_start, overlap_end);
                    // std::cout << "Primary mismatch rate: " << primary_mismatch_rate << std::endl;
                    // std::cout << "Supplementary mismatch rate: " << supp_mismatch_rate << std::endl;

                    // Trim the overlap from the alignment with the higher mismatch
                    // rate
                    if (primary_mismatch_rate > supp_mismatch_rate) {

                        // Handle each scenario
                        if (best_scenario == PRIM_SUPP_BD || best_scenario == PRIM_SUPP_GAP) {
                            // Primary is first, incorporate the overlap into
                            // the beginning of the deletion
                            sv_start -= overlap_length;
                        } else if (best_scenario == SUPP_PRIM_BD || best_scenario == SUPP_PRIM_GAP) {
                            // Primary is last, incorporate the overlap into
                            // the end of the deletion
                            sv_end += overlap_length;
                        }
                    } else {

                        // Handle each scenario
                        if (best_scenario == SUPP_PRIM_BD || best_scenario == SUPP_PRIM_GAP) {
                            // Supplementary is first, incorporate the overlap into
                            // the beginning of the deletion
                            sv_start -= overlap_length;
                        } else if (best_scenario == PRIM_SUPP_BD || best_scenario == PRIM_SUPP_GAP) {
                            // Supplementary is last, incorporate the overlap into
                            // the end of the deletion
                            sv_end += overlap_length;
                        }
                    }
                }
            }

            // Add the best split alignment as the SV call
            // sv_calls.add(primary_chr, sv_start, sv_end, best_supp_type, ".",
            // "SPLITREAD", "./.", best_split_aln_lh_norm);
            std::string sv_type_str = getSVTypeString(best_supp_type);
            sv_count++;
        } else {
            // Resolve complex SVs

            // Simplest case: Largest supplementary is also the closest
            if (largest_supp_alignment == closest_supp_alignment) {
                // [primary] -- [supp_start] -- [supp_end]
                // Determine if opposite strands
                bool opposite_strands = std::get<7>(largest_supp_alignment) != std::get<7>(primary_alignment);

                // Determine if the supplementary alignment is an inversion
                if (opposite_strands) {
                    if (largest_supp_type == SVType::NEUTRAL) {
                        largest_supp_type = SVType::INV;
                    } else if (largest_supp_type == SVType::DUP) {
                        largest_supp_type = SVType::INV_DUP;
                    }
                }

                // Get the SV type strings
                std::string primary_type_str = getSVTypeString(primary_type);
                std::string supp_type_str = getSVTypeString(largest_supp_type);

                // Determine the order of the primary and supplementary
                // alignment to resolve the SV
                if (std::get<1>(largest_supp_alignment) < primary_start) {
                    // [supp_start] -- [supp_end] -- [primary]
                    std::string complex_sv_type_str = supp_type_str + "+" + primary_type_str;

                    // Add the complex SV call
                    addSVCall(sv_calls, (uint32_t)std::get<1>(largest_supp_alignment), (uint32_t)primary_end, "COMPLEX", ".", complex_sv_type_str, "./.", complex_lh_norm);
                    // sv_calls.insert(SVCall{(uint32_t)std::get<1>(largest_supp_alignment), (uint32_t)primary_end, "COMPLEX", ".", complex_sv_type_str, "./.", complex_lh_norm});
                    // sv_calls.add(primary_chr, std::get<1>(largest_supp_alignment), primary_end, SVType::COMPLEX, ".", complex_sv_type_str, "./.", complex_lh_norm);
                    sv_count++;
                } else {
                    // [primary] -- [supp_start] -- [supp_end]
                    std::string complex_sv_type_str = primary_type_str + "+" + supp_type_str;

                    // Add the complex SV call
                    addSVCall(sv_calls, (uint32_t)primary_start, (uint32_t)std::get<2>(largest_supp_alignment), "COMPLEX", ".", complex_sv_type_str, "./.", complex_lh_norm);
                    // sv_calls.insert(SVCall{(uint32_t)primary_start, (uint32_t)std::get<2>(largest_supp_alignment), "COMPLEX", ".", complex_sv_type_str, "./.", complex_lh_norm});
                    // sv_calls.add(primary_chr, primary_start, std::get<2>(largest_supp_alignment), SVType::COMPLEX, ".", complex_sv_type_str, "./.", complex_lh_norm);
                    sv_count++;
                }
            } else {
                // Resolve complex SVs with multiple supplementary alignments
                // Determine the order of the primary and supplementary
                // alignments
                // [primary] -- [closest_supp] -- [largest_supp]
                // [closest_supp] -- [primary] -- [largest_supp]
                // [largest_supp] -- [closest_supp] -- [primary]
                // [largest_supp] -- [primary] -- [closest_supp]
                // Only consider case 1 for efficiency:
                if (primary_end < std::get<1>(closest_supp_alignment) && std::get<2>(closest_supp_alignment) < std::get<1>(largest_supp_alignment)) {
                    // [primary] -- [closest_supp] -- [largest_supp]
                    // Determine if the closest supplementary alignment is an
                    // inversion
                    if (std::get<7>(closest_supp_alignment) != std::get<7>(primary_alignment)) {
                        if (closest_supp_type == SVType::NEUTRAL) {
                            closest_supp_type = SVType::INV;
                        } else if (closest_supp_type == SVType::DUP) {
                            closest_supp_type = SVType::INV_DUP;
                        }
                    }

                    // Run copy number variant predictions on the region between
                    // the closest supplementary alignment and the largest
                    // supplementary alignment
                    SVCandidate sv_candidate(std::get<2>(closest_supp_alignment)+1, std::get<1>(largest_supp_alignment)+1, ".");
                    // std::cout << "TEST4" << std::endl;
                    std::tuple<double, SVType, std::string, bool> result = cnv_caller.runCopyNumberPrediction(primary_chr, sv_candidate);
                    // double complex_log_likelihood = std::get<0>(result);
                    SVType complex_type = std::get<1>(result);

                    // if (std::get<7>(largest_supp_alignment) != std::get<7>(primary_alignment)) {
                    //     if (largest_supp_type == SVType::NEUTRAL) {
                    //         largest_supp_type = SVType::INV;
                    //     } else if (largest_supp_type == SVType::DUP) {
                    //         largest_supp_type = SVType::INV_DUP;
                    //     }
                    // }

                    std::string primary_type_str = getSVTypeString(primary_type);
                    std::string closest_supp_type_str = getSVTypeString(closest_supp_type);
                    // std::string largest_supp_type_str = getSVTypeString(largest_supp_type);
                    // std::string complex_sv_type_str = primary_type_str + "+" + closest_supp_type_str;


                    // Combine the types if equal and not unknown/neutral
                    std::cout << "Resolving complex SVs..." << std::endl;
                    std::string complex_sv_type_str = "";
                    if (primary_type != SVType::UNKNOWN && primary_type != SVType::NEUTRAL) {
                        complex_sv_type_str += primary_type_str;
                        std::cout << "[1] Updated to type: " << complex_sv_type_str << std::endl;
                    }
                    if (closest_supp_type != primary_type && closest_supp_type != SVType::UNKNOWN && closest_supp_type != SVType::NEUTRAL) {
                        if (complex_sv_type_str != "") {
                            complex_sv_type_str += "+";
                        }
                        complex_sv_type_str += closest_supp_type_str;
                        std::cout << "[2] Updated to type: " << complex_sv_type_str << std::endl;
                    }
                    if (complex_type != closest_supp_type && complex_type != primary_type && complex_type != SVType::UNKNOWN && complex_type != SVType::NEUTRAL) {
                        if (complex_sv_type_str != "") {
                            complex_sv_type_str += "+";
                        }
                        complex_sv_type_str += getSVTypeString(complex_type);
                        std::cout << "[3] Updated to type: " << complex_sv_type_str << std::endl;
                    }

                    // Add the complex SV call if not empty
                    if (complex_sv_type_str != "") {
                        std::cout << "Found complex SV type: " << complex_sv_type_str << std::endl;
                        // sv_calls.add(primary_chr, primary_start,
                        // std::get<2>(largest_supp_alignment), SVType::COMPLEX,
                        // ".", complex_sv_type_str, "./.", complex_lh_norm);
                        // sv_calls.insert(SVCall{(uint32_t)primary_start, (uint32_t)std::get<2>(largest_supp_alignment), "COMPLEX", ".", complex_sv_type_str, "./.", complex_lh_norm});
                        addSVCall(sv_calls, (uint32_t)primary_start, (uint32_t)std::get<2>(largest_supp_alignment), "COMPLEX", ".", complex_sv_type_str, "./.", complex_lh_norm);
                        sv_count++;
                    }
                }                
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
    // std::set<std::string> chrs = this->getChromosomes();
   //for (auto const& chr : chrs) {
   for (const auto& pair : sv_calls) {
        // if (this->sv_calls.find(chr) == this->sv_calls.end()) {
        //     continue;
        // }
        std::string chr = pair.first;
        const std::set<SVCall>& sv_calls = pair.second;
        std::cout << "Saving SV calls for " << chr << "..." << std::endl;
        // for (auto const& sv_call : this->sv_calls[chr]) {
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

