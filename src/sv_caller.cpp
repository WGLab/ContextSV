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

#include "utils.h"
#include "sv_types.h"
/// @endcond

# define DUP_SEQSIM_THRESHOLD 0.9  // Sequence similarity threshold for duplication detection

int SVCaller::readNextAlignment(samFile *fp_in, hts_itr_t *itr, bam1_t *bam1)
{
    int ret = sam_itr_next(fp_in, itr, bam1);
    return ret;
}

RegionData SVCaller::detectSVsFromRegion(std::string region)
{
    // Open the BAM file
    std::string bam_filepath = this->input_data->getLongReadBam();
    samFile *fp_in = sam_open(bam_filepath.c_str(), "r");
    if (fp_in == NULL) {
        std::cerr << "ERROR: failed to open " << bam_filepath << std::endl;
        exit(1);
    }

    // Load the header for the BAM file
    bam_hdr_t *bamHdr = sam_hdr_read(fp_in);
    if (bamHdr == NULL) {
        std::cerr << "ERROR: failed to read header for " << bam_filepath << std::endl;
        exit(1);
    }

    // Load the index for the BAM file
    hts_idx_t *idx = sam_index_load(fp_in, bam_filepath.c_str());
    if (idx == NULL) {
        std::cerr << "ERROR: failed to load index for " << bam_filepath << std::endl;
        exit(1);
    }

    // Create a read and iterator for the region
    bam1_t *bam1 = bam_init1();
    hts_itr_t *itr = sam_itr_querys(idx, bamHdr, region.c_str());

    // Main loop to process the alignments
    SVData sv_calls;
    int num_alignments = 0;
    PrimaryMap primary_alignments;
    SuppMap supplementary_alignments;
    while (readNextAlignment(fp_in, itr, bam1) >= 0) {

        // Skip secondary and unmapped alignments, duplicates, and QC failures
        if (bam1->core.flag & BAM_FSECONDARY || bam1->core.flag & BAM_FUNMAP || bam1->core.flag & BAM_FDUP || bam1->core.flag & BAM_FQCFAIL) {
            // Do nothing

        // Skip alignments with low mapping quality
        } else if (bam1->core.qual < this->min_mapq) {
            // Do nothing

        } else {
            std::string qname = bam_get_qname(bam1);  // Query template name

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
                AlignmentData alignment(chr, start, end, ".", query_start, query_end, std::move(match_map), fwd_strand);
                primary_alignments[qname] = std::move(alignment);

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
                AlignmentData alignment(chr, start, end, ".", query_start, query_end, std::move(match_map), fwd_strand);
                supplementary_alignments[qname].emplace_back(alignment);
            }
        }

        num_alignments++;
    }

    hts_itr_destroy(itr);
    bam_destroy1(bam1);
    sam_close(fp_in);
    bam_hdr_destroy(bamHdr);
    hts_idx_destroy(idx);

    // Return the SV calls and the primary and supplementary alignments
    return std::make_tuple(std::move(sv_calls), std::move(primary_alignments), std::move(supplementary_alignments));
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

SVCaller::SVCaller(InputData &input_data)
{
    this->input_data = &input_data;
}

std::tuple<std::unordered_map<int, int>, int32_t, int32_t> SVCaller::detectSVsFromCIGAR(bam_hdr_t* header, bam1_t* alignment, SVData& sv_calls, bool is_primary)
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
                    std::string window_str = this->input_data->queryRefGenome(chr, ins_ref_pos, ins_ref_pos + op_len - 1);

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

                // Add to SV calls (1-based) with the appropriate SV type
                ref_pos = pos+1;
                ref_end = ref_pos + op_len -1;
                if (is_duplication) {
                    sv_calls.add(chr, ref_pos, ref_end, SVType::DUP, ins_seq_str, "CIGARDUP", "./.", 0.0);
                } else {
                    sv_calls.add(chr, ref_pos, ref_end, SVType::INS, ins_seq_str, "CIGARINS", "./.", 0.0);
                }
            }

        // Check if the CIGAR operation is a deletion
        } else if (op == BAM_CDEL && is_primary) {

            // Add the SV if greater than the minimum SV size
            if (op_len >= this->min_sv_size)
            {
                ref_pos = pos+1;
                ref_end = ref_pos + op_len -1;
                sv_calls.add(chr, ref_pos, ref_end, SVType::DEL, ".", "CIGARDEL", "./.", 0.0);  // Add the deletion
            }

        // Check if the CIGAR operation is a clipped base
        } else if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) {

            sv_calls.updateClippedBaseSupport(chr, pos);  // Update clipped base support

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
            std::string cmatch_ref_str = this->input_data->queryRefGenome(chr, cmatch_pos, cmatch_pos + op_len - 1);

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

SVData SVCaller::run()
{
    // Get the chromosomes to process
    std::vector<std::string> chromosomes;
    if (this->input_data->getChromosome() != "") {
        chromosomes.push_back(this->input_data->getChromosome());
    } else {
        chromosomes = this->input_data->getRefGenomeChromosomes();
    }
    int chr_count = chromosomes.size();

    // Loop through each region and detect SVs in chunks
    std::cout << "Detecting SVs from " << chr_count << " chromosome(s)..." << std::endl;
    int chunk_count = 100;  // Number of chunks to split the chromosome into
    int region_count = 0;
    SVData sv_calls;
    int min_cnv_length = this->input_data->getMinCNVLength();
    for (const auto& chr : chromosomes) {
        std::cout << "Running SV detection for chromosome " << chr << "..." << std::endl;

        // Split the chromosome into chunks
        std::vector<std::string> region_chunks;
        if (this->input_data->isRegionSet()) {

            // Use one chunk for the specified region
            std::pair<int32_t, int32_t> region = this->input_data->getRegion();
            int region_start = region.first;
            int region_end = region.second;
            std::string chunk = chr + ":" + std::to_string(region_start) + "-" + std::to_string(region_end);
            region_chunks.push_back(chunk);
            std::cout << "Using specified region " << chunk << "..." << std::endl;
            
        } else {
            int chr_len = this->input_data->getRefGenomeChromosomeLength(chr);
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
        CNVCaller cnv_caller(*this->input_data);
        cnv_caller.loadChromosomeData(chr);

        // Process each chunk one at a time
        std::cout << "Processing " << region_chunks.size() << " region(s) for chromosome " << chr << "..." << std::endl;
        for (const auto& sub_region : region_chunks) {
            // std::cout << "Detecting CIGAR string SVs from " << sub_region << "..." << std::endl;
            RegionData region_data = this->detectSVsFromRegion(sub_region);
            SVData& sv_calls_region = std::get<0>(region_data);
            PrimaryMap& primary_map = std::get<1>(region_data);
            SuppMap& supp_map = std::get<2>(region_data);
            int region_sv_count = sv_calls_region.totalCalls();
            if (region_sv_count > 0) {
                std::cout << "Detected " << region_sv_count << " SVs from " << sub_region << "..." << std::endl;
            }

            // Run copy number variant predictions on the SVs detected from the
            // CIGAR string, using a minimum CNV length threshold
            // std::cout << "Detecting copy number variants from CIGAR string SVs..." << std::endl;
            std::map<SVCandidate, SVInfo>& cigar_svs = sv_calls_region.getChromosomeSVs(chr);
            if (cigar_svs.size() > 0) {
                std::cout << "Running copy number variant detection from CIGAR string SVs..." << std::endl;
                cnv_caller.runCIGARCopyNumberPrediction(chr, cigar_svs, min_cnv_length);
            }

            // Run split-read SV detection in a single thread, combined with
            // copy number variant predictions
            std::cout << "Detecting copy number variants from split reads..." << std::endl;
            this->detectSVsFromSplitReads(sv_calls_region, primary_map, supp_map, cnv_caller);
            sv_calls.concatenate(sv_calls_region);  // Add the calls to the main set
        }

        region_count++;
        std::cout << "Completed " << region_count << " of " << chr_count << " chromosome(s)..." << std::endl;
    }

    std::cout << "SV calling completed." << std::endl;

    return sv_calls;
}


// Detect SVs from split read alignments
void SVCaller::detectSVsFromSplitReads(SVData& sv_calls, PrimaryMap& primary_map, SuppMap& supp_map, CNVCaller& cnv_caller)
{
    // Find split-read SV evidence
    int sv_count = 0;
    int min_cnv_length = this->input_data->getMinCNVLength();
    for (const auto& entry : primary_map) {
        std::string qname = entry.first;
        AlignmentData primary_alignment = entry.second;
        std::string primary_chr = std::get<0>(primary_alignment);
        int32_t primary_start = std::get<1>(primary_alignment);
        int32_t primary_end = std::get<2>(primary_alignment);
        int32_t primary_query_start = std::get<4>(primary_alignment);
        int32_t primary_query_end = std::get<5>(primary_alignment);
        std::unordered_map<int, int> primary_match_map = std::get<6>(primary_alignment);
        bool primary_strand = std::get<7>(primary_alignment);
        if (supp_map.find(qname) == supp_map.end()) {
            continue;
        }

        // Resolve overlaps between the primary and supplementary query
        // sequences
        for (auto it = supp_map[qname].begin(); it != supp_map[qname].end(); ++it) {
            std::string supp_chr = std::get<0>(*it);
            // int32_t supp_start = std::get<1>(*it);
            // int32_t supp_end = std::get<2>(*it);
            int32_t supp_query_start = std::get<4>(*it);
            int32_t supp_query_end = std::get<5>(*it);
            std::unordered_map<int, int> supp_match_map = std::get<6>(*it);
            bool supp_strand = std::get<7>(*it);

            // Resolve overlaps between the primary and supplementary query
            // sequences
            int32_t overlap_start = std::max(primary_query_start, supp_query_start);
            int32_t overlap_end = std::min(primary_query_end, supp_query_end);
            int32_t overlap_length = overlap_end - overlap_start;
            if (overlap_length > 0) {

                // Calculate the mismatch rate for each alignment at the overlap
                double primary_mismatch_rate = this->calculateMismatchRate(primary_match_map, overlap_start, overlap_end-1);
                double supp_mismatch_rate = this->calculateMismatchRate(supp_match_map, overlap_start, overlap_end-1);
                // std::cout << "Primary mismatch rate: " << primary_mismatch_rate << std::endl;
                // std::cout << "Supplementary mismatch rate: " << supp_mismatch_rate << std::endl;

                // Trim the overlap from the alignment with the higher mismatch
                // rate
                if (primary_mismatch_rate > supp_mismatch_rate) {
                    if (overlap_start == primary_query_start) {
                        primary_start += overlap_length;
                    } else if (overlap_end == primary_query_end) {
                        primary_end -= overlap_length;
                    }

                } else {
                    if (overlap_start == supp_query_start) {
                        // supp_start += overlap_length;
                        // Update the value in the supp map
                        std::get<1>(*it) += overlap_length;
                    } else if (overlap_end == supp_query_end) {
                        // supp_end -= overlap_length;
                        // Update the value in the supp map
                        std::get<2>(*it) -= overlap_length;
                    }
                }
            }
        }

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
        double primary_log_likelihood = std::numeric_limits<double>::lowest();
        if (primary_end - primary_start >= min_cnv_length) {
            SVCandidate sv_candidate(primary_start+1, primary_end+1, ".");
            std::tuple<double, SVType, std::string, bool> result = cnv_caller.runCopyNumberPrediction(primary_chr, sv_candidate);
            primary_log_likelihood = std::get<0>(result);
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
            const auto& supp_chr = std::get<0>(*it);
            int32_t supp_start = std::get<1>(*it);
            int32_t supp_end = std::get<2>(*it);
            int32_t supp_length = supp_end - supp_start;
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
        double largest_supp_log_likelihood = std::numeric_limits<double>::lowest();
        SVType largest_supp_type = SVType::UNKNOWN;
        if (largest_supp_length >= min_cnv_length) {
            SVCandidate sv_candidate(std::get<1>(largest_supp_alignment)+1, std::get<2>(largest_supp_alignment)+1, ".");
            std::tuple<double, SVType, std::string, bool> result = cnv_caller.runCopyNumberPrediction(primary_chr, sv_candidate);
            largest_supp_log_likelihood = std::get<0>(result);
            // largest_supp_log_likelihood /= (double)largest_supp_length;  // Normalize the log likelihood by the length
            largest_supp_type = std::get<1>(result);
        }

        // Run copy number variant predictions on the closest non-overlapping
        // supplementary alignment (if not the same as the largest)
        double closest_supp_log_likelihood = std::numeric_limits<double>::lowest();
        SVType closest_supp_type = SVType::UNKNOWN;
        if (largest_supp_alignment != closest_supp_alignment) {
            if (closest_supp_length >= min_cnv_length) {
                SVCandidate sv_candidate(std::get<1>(closest_supp_alignment)+1, std::get<2>(closest_supp_alignment)+1, ".");
                std::tuple<double, SVType, std::string, bool> result = cnv_caller.runCopyNumberPrediction(primary_chr, sv_candidate);
                closest_supp_log_likelihood = std::get<0>(result);
                // closest_supp_log_likelihood /= (double)closest_supp_length;  // Normalize the log likelihood by the length
                closest_supp_type = std::get<1>(result);
            }
        }

        // Loop through all the supplementary alignments and find the highest
        // likelihood prediction
        double best_supp_log_likelihood = std::numeric_limits<double>::lowest();
        SVType best_supp_type = SVType::UNKNOWN;
        std::pair<int32_t, int32_t> best_supp_candidate;
        for (auto it = supp_map[qname].begin(); it != supp_map[qname].end(); ++it) {
            int32_t supp_start = std::get<1>(*it);
            int32_t supp_end = std::get<2>(*it);

            // Create the SV candidate as the boundary of the supplementary
            // and primary alignment
            // int32_t sv_start = std::min(primary_start, std::get<1>(*it));
            // int32_t sv_end = std::max(primary_end, std::get<2>(*it));
            int32_t sv_start = std::min(primary_start, supp_start);
            int32_t sv_end = std::max(primary_end, supp_end);
            SVCandidate sv_candidate(sv_start+1, sv_end+1, ".");

            // Determine if the strand is the same as the primary alignment
            bool same_strand = std::get<7>(*it) == std::get<7>(primary_alignment);

            // SVCandidate sv_candidate(std::get<1>(*it)+1, std::get<2>(*it)+1, ".");
            std::tuple<double, SVType, std::string, bool> result = cnv_caller.runCopyNumberPrediction(primary_chr, sv_candidate);
            double supp_likelihood = std::get<0>(result);
            SVType supp_type = std::get<1>(result);

            // If opposite strand, set the type to INV or INV_DUP
            if (!same_strand) {
                if (supp_type == SVType::NEUTRAL) {
                    supp_type = SVType::INV;
                } else if (supp_type == SVType::DUP) {
                    supp_type = SVType::INV_DUP;
                }
            }

            if (supp_type != SVType::UNKNOWN && supp_likelihood > best_supp_log_likelihood) {
                best_supp_log_likelihood = supp_likelihood;
                // best_supp_log_likelihood /= (double)(sv_end - sv_start);  // Normalize the log likelihood by the length
                best_supp_type = supp_type;
                best_supp_candidate = std::make_pair(supp_start, supp_end);
            }
        }

        // Add the SV call with the highest likelihood prediction
        if (best_supp_log_likelihood > primary_log_likelihood || best_supp_log_likelihood > largest_supp_log_likelihood || best_supp_log_likelihood > closest_supp_log_likelihood) {
            int32_t sv_start = best_supp_candidate.first;
            int32_t sv_end = best_supp_candidate.second;
            sv_calls.add(primary_chr, sv_start, sv_end, best_supp_type, ".", "SPLITREAD", "./.", best_supp_log_likelihood);
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
                    sv_calls.add(primary_chr, std::get<1>(largest_supp_alignment), primary_end, SVType::COMPLEX, ".", complex_sv_type_str, "./.", 0.0);
                    sv_count++;
                } else {
                    // [primary] -- [supp_start] -- [supp_end]
                    std::string complex_sv_type_str = primary_type_str + "+" + supp_type_str;

                    // Add the complex SV call
                    sv_calls.add(primary_chr, primary_start, std::get<2>(largest_supp_alignment), SVType::COMPLEX, ".", complex_sv_type_str, "./.", 0.0);
                    sv_count++;
                }
            }
        }
    }

    // Print the number of SVs detected from split-read alignments
    if (sv_count > 0) {
        std::cout << "Found " << sv_count << " SVs from split-read alignments" << std::endl;
    }
}
