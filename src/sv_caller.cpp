//
// sv_caller.cpp:
// Detect SVs from long read alignments
//

#include "sv_caller.h"
#include "cnv_caller.h"

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

#include "utils.h"
#include "sv_types.h"
/// @endcond

# define DUP_SEQSIM_THRESHOLD 0.9  // Sequence similarity threshold for duplication detection

int SVCaller::readNextAlignment(samFile *fp_in, hts_itr_t *itr, bam1_t *bam1)
{
    // Read the next alignment
    int ret = sam_itr_next(fp_in, itr, bam1);

    // Return the result of reading the next alignment
    return ret;
}

RegionData SVCaller::detectSVsFromRegion(std::string region)
{
    SVData sv_calls;
    std::string bam_filepath = this->input_data->getLongReadBam();

    // Open the BAM file in a thread-safe manner
    samFile *fp_in = sam_open(bam_filepath.c_str(), "r");
    if (fp_in == NULL) {
        std::cerr << "ERROR: failed to open " << bam_filepath << std::endl;
        exit(1);
    }

    // Get the header in a thread-safe manner
    bam_hdr_t *bamHdr = sam_hdr_read(fp_in);
    if (bamHdr == NULL) {
        std::cerr << "ERROR: failed to read header for " << bam_filepath << std::endl;
        exit(1);
    }

    // Get the index in a thread-safe manner
    hts_idx_t *idx = sam_index_load(fp_in, bam_filepath.c_str());
    if (idx == NULL) {
        std::cerr << "ERROR: failed to load index for " << bam_filepath << std::endl;
        exit(1);
    }

    // Create a read and iterator for the region in a thread-safe manner
    bam1_t *bam1 = bam_init1();
    hts_itr_t *itr = sam_itr_querys(idx, bamHdr, region.c_str());

    // Loop through the alignments
    // Create a map of primary and supplementary alignments by QNAME (query template name)
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
            // Get the QNAME (query template name) for associating split reads
            std::string qname = bam_get_qname(bam1);

            // Process primary alignments
            if (!(bam1->core.flag & BAM_FSUPPLEMENTARY)) {

                // Get the primary alignment chromosome, start, end, and depth
                std::string chr = bamHdr->target_name[bam1->core.tid];
                int64_t start = bam1->core.pos;
                int64_t end = bam_endpos(bam1);

                // Call SVs directly from the CIGAR string
                std::tuple<std::unordered_map<int, int>, int32_t, int32_t> query_info = this->detectSVsFromCIGAR(bamHdr, bam1, sv_calls, true);
                std::unordered_map<int, int> match_map = std::get<0>(query_info);
                int32_t query_start = std::get<1>(query_info);
                int32_t query_end = std::get<2>(query_info);

                // Add the primary alignment to the map
                AlignmentData alignment(chr, start, end, ".", query_start, query_end, match_map);
                primary_alignments[qname] = std::move(alignment);

            // Process supplementary alignments
            } else if (bam1->core.flag & BAM_FSUPPLEMENTARY) {

                // Add the supplementary alignment to the map
                std::string chr = bamHdr->target_name[bam1->core.tid];
                int32_t start = bam1->core.pos;
                int32_t end = bam_endpos(bam1);

                // Get CIGAR string information, but don't call SVs
                std::tuple<std::unordered_map<int, int>, int32_t, int32_t> query_info = this->detectSVsFromCIGAR(bamHdr, bam1, sv_calls, false);
                const std::unordered_map<int, int>& match_map = std::get<0>(query_info);
                int32_t query_start = std::get<1>(query_info);
                int32_t query_end = std::get<2>(query_info);

                // Add the supplementary alignment to the map
                AlignmentData alignment(chr, start, end, ".", query_start, query_end, std::move(match_map));
                supplementary_alignments[qname].emplace_back(alignment);
            }
        }

        // Increment the number of alignment records processed
        num_alignments++;
    }

    // Destroy the iterator
    hts_itr_destroy(itr);

    // Destroy the read
    bam_destroy1(bam1);

    // Close the BAM file
    sam_close(fp_in);

    // Destroy the header
    bam_hdr_destroy(bamHdr);

    // Destroy the index
    hts_idx_destroy(idx);

    // Return the SV calls and the primary and supplementary alignments
    // return std::make_tuple(sv_calls, primary_alignments,
    // supplementary_alignments);
    return std::make_tuple(std::move(sv_calls), std::move(primary_alignments), std::move(supplementary_alignments));
}

double SVCaller::calculateMismatchRate(std::unordered_map<int, int> &match_map, int32_t start, int32_t end)
{
    // Calculate the mismatch rate
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

    // Calculate the mismatch rate
    double mismatch_rate = (double)mismatch_count / (double)(match_count + mismatch_count);

    // Return the mismatch rate
    return mismatch_rate;
}

SVCaller::SVCaller(InputData &input_data)
{
    this->input_data = &input_data;
}

std::tuple<std::unordered_map<int, int>, int32_t, int32_t> SVCaller::detectSVsFromCIGAR(bam_hdr_t* header, bam1_t* alignment, SVData& sv_calls, bool is_primary)
{
    // Get the chromosome
    std::string chr = header->target_name[alignment->core.tid];

    // Get the position of the alignment in the reference genome
    int32_t pos = alignment->core.pos;

    // Get the CIGAR string
    uint32_t* cigar = bam_get_cigar(alignment);

    // Get the CIGAR length
    int cigar_len = alignment->core.n_cigar;

    // Track the query position
    int query_pos = 0;

    // Loop through the CIGAR string (0-based) and detect insertions and deletions in
    // reference coordinates (1-based)
    // POS is the leftmost position of where the alignment maps to the reference:
    // https://genome.sph.umich.edu/wiki/SAM
    // std::vector<std::thread> threads;
    // std::vector<SVData> sv_calls_vec;

    // Create a map of query position to match/mismatch (1/0) for calculating
    // the mismatch rate at alignment overlaps
    std::unordered_map<int, int> query_match_map;

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

        // Get the CIGAR operation
        int op = bam_cigar_op(cigar[i]);

        // Get the CIGAR operation length
        int op_len = bam_cigar_oplen(cigar[i]);
        
        // Check if the CIGAR operation is an insertion
        if (op == BAM_CINS && is_primary) {

            // Add the SV if greater than the minimum SV size
            if (op_len >= this->min_sv_size) {

                // Get the sequence of the insertion from the query
                // std::string ins_seq_str = "";
                // uint8_t* seq_ptr = bam_get_seq(alignment);
                // for (int j = 0; j < op_len; j++) {
                //     ins_seq_str += seq_nt16_str[bam_seqi(seq_ptr, query_pos + j)];
                // }
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

                // Lock the SV calls object and add the insertion
                std::lock_guard<std::mutex> lock(this->sv_mtx);
                if (is_duplication) {
                    sv_calls.add(chr, ref_pos, ref_end, DUP, ins_seq_str, "CIGARDUP");
                } else {
                    sv_calls.add(chr, ref_pos, ref_end, INS, ins_seq_str, "CIGARINS");
                }
            }

        // Check if the CIGAR operation is a deletion
        } else if (op == BAM_CDEL && is_primary) {

            // Add the SV if greater than the minimum SV size
            if (op_len >= this->min_sv_size) {
                
                // Add the deletion to the SV calls (1-based)
                ref_pos = pos+1;
                ref_end = ref_pos + op_len -1;

                // Lock the SV calls object and add the deletion
                // std::lock_guard<std::mutex> lock(this->sv_mtx);
                sv_calls.add(chr, ref_pos, ref_end, DEL, ".", "CIGARDEL");
            }

        // Check if the CIGAR operation is a soft clip
        } else if (op == BAM_CSOFT_CLIP) {

            // Update the clipped base support
            // std::lock_guard<std::mutex> lock(this->sv_mtx);
            sv_calls.updateClippedBaseSupport(chr, pos);

            // Update the query alignment start position
            if (!first_op) {
                query_start = query_pos + op_len;
                first_op = true;
            }
        }

        // Update match/mismatch query map
        if (op == BAM_CEQUAL) {
            // match_count += op_len;
            for (int j = 0; j < op_len; j++) {
                query_match_map[query_pos + j] = 1;
            }
        } else if (op == BAM_CDIFF) {
            // mismatch_count += op_len;
            for (int j = 0; j < op_len; j++) {
                query_match_map[query_pos + j] = 0;
            }
        } else if (op == BAM_CMATCH) {
            // Compare read and reference sequences
            // Get the sequence from the query
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

    // Update the query end position
    query_end = query_pos;

    // Return the mismatch map and the query start and end positions
    return std::tuple<std::unordered_map<int, int>, int32_t, int32_t>(query_match_map, query_start, query_end);
}

// Detect SVs from split read alignments (primary and supplementary) and
// directly from the CIGAR string
SVData SVCaller::run()
{
    // Open the BAM file
    std::string bam_filepath = this->input_data->getLongReadBam();

    // Get the region data
    std::vector<std::string> chromosomes;
    if (this->input_data->getChromosome() != "") {
        chromosomes.push_back(this->input_data->getChromosome());
    } else {
        chromosomes = this->input_data->getRefGenomeChromosomes();
    }
    int chr_count = chromosomes.size();

    // Loop through each region and detect SVs
    std::cout << "Detecting SVs from " << chr_count << " chromosome(s)..." << std::endl;
    int region_count = 0;
    auto start1 = std::chrono::high_resolution_clock::now();
    SVData sv_calls;
    int chunk_count = 10000;  // Number of chunks to split the chromosome into
    for (const auto& chr : chromosomes) {
        std::cout << "Extracting alignments for chromosome " << chr << "..." << std::endl;

        // Split the chromosome into chunks
        std::vector<std::string> region_chunks;

        // Get the region start and end positions
        if (this->input_data->isRegionSet()) {
            std::pair<int32_t, int32_t> region = this->input_data->getRegion();
            int region_start = region.first;
            int region_end = region.second;

            // Use one chunk for the region
            std::string chunk = chr + ":" + std::to_string(region_start) + "-" + std::to_string(region_end);
            region_chunks.push_back(chunk);
            
        } else {
            int chr_len = this->input_data->getRefGenomeChromosomeLength(chr);
            int chunk_size = chr_len / chunk_count;
            for (int i = 0; i < chunk_count; i++) {
                int start = i * chunk_size + 1;  // 1-based
                int end = start + chunk_size;
                if (i == chunk_count - 1) {
                    end = chr_len;
                }
                std::string chunk = chr + ":" + std::to_string(start) + "-" + std::to_string(end);
                region_chunks.push_back(chunk);
            }
        }

        // Load chromosome data for copy number predictions
        std::cout << "Loading chromosome data for copy number predictions..." << std::endl;
        CNVCaller cnv_caller(*this->input_data);
        cnv_caller.loadChromosomeData(chr);
        std::cout << "Loaded chromosome data for copy number predictions." << std::endl;

        // Process each chunk one at a time
        for (const auto& sub_region : region_chunks) {
            // Detect SVs from the sub-region
            std::cout << "Detecting SVs from " << sub_region << "..." << std::endl;
            RegionData region_data = this->detectSVsFromRegion(sub_region);
            SVData& sv_calls_region = std::get<0>(region_data);
            PrimaryMap& primary_map = std::get<1>(region_data);
            SuppMap& supp_map = std::get<2>(region_data);
            int region_sv_count = sv_calls_region.totalCalls();
            std::cout << "Detected " << region_sv_count << " SVs from " << sub_region << "..." << std::endl;

            // Run split-read SV detection in a single thread
            std::cout << "Detecting SVs from split-read alignments..." << std::endl;
            this->detectSVsFromSplitReads(sv_calls_region, primary_map, supp_map);
            int split_read_sv_count = sv_calls_region.totalCalls() - region_sv_count;
            std::cout << "Detected " << split_read_sv_count << " additional SVs from split-read alignments..." << std::endl;
            // sv_calls.concatenate(sv_calls_region);

            // Run copy number predictions
            std::cout << "Predicting copy number from split-read alignments..." << std::endl;
            std::map<SVCandidate, SVInfo>& sv_candidates = sv_calls_region.getChromosomeSVs(chr);
            cnv_caller.runCopyNumberPrediction(chr, sv_candidates);
            std::cout << "Predicted copy number from split-read alignments." << std::endl;

            // Add the SV calls to the main SV calls object
            sv_calls.concatenate(sv_calls_region);
        }


        // Increment the region count
        region_count++;
        std::cout << "Extracted aligments for " << region_count << " of " << chr_count << " chromosome(s)..." << std::endl;
    }

    auto end1 = std::chrono::high_resolution_clock::now();
    std::cout << "Finished detecting " << sv_calls.totalCalls() << " SVs from " << chr_count << " chromosome(s). Elapsed time: " << getElapsedTime(start1, end1) << std::endl;

    return sv_calls;
}


// Detect SVs from split read alignments
void SVCaller::detectSVsFromSplitReads(SVData& sv_calls, PrimaryMap& primary_map, SuppMap& supp_map)
{
    // Loop through the map of primary alignments by QNAME and find gaps and
    // overlaps from supplementary alignments
    int sv_count = 0;
    for (const auto& entry : primary_map) {
        // Get the QNAME
        std::string qname = entry.first;

        // Get the first primary alignment
        AlignmentData primary_alignment = entry.second;

        // Get the primary alignment chromosome
        std::string primary_chr = std::get<0>(primary_alignment);

        // Get the start and end positions of the primary alignment (reference)
        int32_t primary_start = std::get<1>(primary_alignment);
        int32_t primary_end = std::get<2>(primary_alignment);

        // Get the query start and end positions
        int32_t primary_query_start = std::get<4>(primary_alignment);
        int32_t primary_query_end = std::get<5>(primary_alignment);

        // Get the primary mismatch map
        std::unordered_map<int, int> primary_match_map = std::get<6>(primary_alignment);

        // Loop through the supplementary alignments and find gaps and overlaps
        AlignmentVector supp_alignments = supp_map[qname];
        for (const auto& supp_alignment : supp_alignments) {

            // Get the supplementary alignment chromosome
            std::string supp_chr = std::get<0>(supp_alignment);

            // Skip supplementary alignments that are on a different chromosome
            // for now (TODO: Use for identifying trans-chromosomal SVs such as
            // translocations)
            if (primary_chr != supp_chr) {
                continue;
            }

            // Get the start and end positions of the supplementary alignment
            // (on the reference genome)
            int32_t supp_start = std::get<1>(supp_alignment);
            int32_t supp_end = std::get<2>(supp_alignment);

            // Get the query start and end positions
            int32_t supp_query_start = std::get<4>(supp_alignment);
            int32_t supp_query_end = std::get<5>(supp_alignment);

            // Get the supplementary mismatch map
            std::unordered_map<int, int> supp_match_map = std::get<6>(supp_alignment);

            // Determine if there is overlap between the primary and
            // supplementary query sequences
            int32_t overlap_start = std::max(primary_query_start, supp_query_start);
            int32_t overlap_end = std::min(primary_query_end, supp_query_end);
            int32_t overlap_length = overlap_end - overlap_start;
            if (overlap_length > 0) {
                // Calculate the mismatch rate at the overlap for each alignment at the overlap
                double primary_mismatch_rate = this->calculateMismatchRate(primary_match_map, overlap_start, overlap_end);
                double supp_mismatch_rate = this->calculateMismatchRate(supp_match_map, overlap_start, overlap_end);

                // Trim the overlap from the alignment with the higher mismatch
                // rate
                if (primary_mismatch_rate > supp_mismatch_rate) {
                    // Trim the overlap from the primary alignment
                    primary_end -= overlap_length;
                } else {
                    // Trim the overlap from the supplementary alignment
                    supp_end -= overlap_length;
                }
            }

            // Gap analysis (deletion or duplication)
            if (supp_start < primary_start && supp_end < primary_start) {
                // Gap with supplementary before primary:
                // [supp_start] [supp_end] -- [primary_start] [primary_end]

                // Use the gap ends as the SV endpoints
                if (primary_start - supp_end >= this->min_sv_size) {

                    // Add the SV call
                    sv_calls.add(supp_chr, supp_end+1, primary_start+1, UNKNOWN, ".", "GAPINNER_A");
                    sv_count++;
                }

                // Also use the alignment ends as the SV endpoints
                if (primary_end - supp_start >= this->min_sv_size) {

                    // Add the SV call
                    sv_calls.add(supp_chr, supp_start+1, primary_end+1, UNKNOWN, ".", "GAPOUTER_A");
                    sv_count++;
                }

                
            } else if (supp_start > primary_end && supp_end > primary_end) {
                // Gap with supplementary after primary:
                // [primary_start] [primary_end] -- [supp_start] [supp_end]

                // Use the gap ends as the SV endpoints
                if (supp_start - primary_end >= this->min_sv_size) {

                    // Add the SV call
                    sv_calls.add(supp_chr, primary_end+1, supp_start+1, UNKNOWN, ".", "GAPINNER_B");
                    sv_count++;
                }

                // Also use the alignment ends as the SV endpoints
                if (supp_end - primary_start >= this->min_sv_size) {

                    // Add the SV call
                    sv_calls.add(supp_chr, primary_start+1, supp_end+1, UNKNOWN, ".", "GAPOUTER_B");
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
