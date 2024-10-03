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
                std::tuple<double, int32_t, int32_t> query_info = this->detectSVsFromCIGAR(bamHdr, bam1, sv_calls, true);
                double mismatch_rate = std::get<0>(query_info);
                int32_t query_start = std::get<1>(query_info);
                int32_t query_end = std::get<2>(query_info);

                // Add the primary alignment to the map
                AlignmentData alignment(chr, start, end, ".", query_start, query_end, mismatch_rate);
                // primary_alignments[qname].push_back(alignment);
                primary_alignments[qname] = alignment;

                // Check if there are supplementary alignments, determine if
                // there is overlap, and trim the alignments with the highest
                // mismatch rate
                auto it = supplementary_alignments.find(qname);
                if (it != supplementary_alignments.end()) {
                    AlignmentVector supp_alignments = it->second;
                    for (auto& supp_alignment : supp_alignments) {
                        // Get the query start and end positions
                        int32_t supp_query_start = std::get<4>(supp_alignment);
                        int32_t supp_query_end = std::get<5>(supp_alignment);

                        // Check if there is overlap between the primary and
                        // supplementary alignments
                        int overlap = std::max(0, std::min(query_end, supp_query_end) - std::max(query_start, supp_query_start));
                        if (overlap > 0) {
                            // Calculate the mismatch rate at the overlap
                            query_ovelap_start = std::max(query_start, supp_query_start);
                            query_overlap_end = std::min(query_end, supp_query_end);

                            // Get the sequence of the primary alignment
                            std::string primary_seq = bam_get_seq(bam1);
                            std::string supp_seq = bam_get_seq(supp_alignment);

                            // TODO: Calculate the mismatch rate at the overlap
                            int mismatch_count = 0;
                            for (int i = query_overlap_start; i < query_overlap_end; i++) {
                                if (primary_seq[i] != supp_seq[i]) {
                                    mismatch_count++;
                                }
                            }
                        }
                    }
                }

            // Process supplementary alignments
            } else if (bam1->core.flag & BAM_FSUPPLEMENTARY) {

                // Add the supplementary alignment to the map
                std::string chr = bamHdr->target_name[bam1->core.tid];
                int32_t start = bam1->core.pos;
                int32_t end = bam_endpos(bam1);

                // Get CIGAR string information (don't call SVs in the supplementary alignments)
                std::tuple<double, int32_t, int32_t> query_info = this->detectSVsFromCIGAR(bamHdr, bam1, sv_calls, false);
                double mismatch_rate = std::get<0>(query_info);
                int32_t query_start = std::get<1>(query_info);
                int32_t query_end = std::get<2>(query_info);

                // Add the supplementary alignment to the map
                AlignmentData alignment(chr, start, end, ".", query_start, query_end, mismatch_rate);
                auto it = supplementary_alignments.find(qname);
                if (it == supplementary_alignments.end()) {
                    supplementary_alignments[qname] = {alignment};
                } else {
                    supplementary_alignments[qname].push_back(alignment);
                }
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
    return std::make_tuple(sv_calls, primary_alignments, supplementary_alignments);
}

SVCaller::SVCaller(InputData &input_data)
{
    this->input_data = &input_data;
}
std::tuple<double, int32_t, int32_t> SVCaller::detectSVsFromCIGAR(bam_hdr_t* header, bam1_t* alignment, SVData& sv_calls, bool is_primary)
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
    std::vector<std::thread> threads;
    std::vector<SVData> sv_calls_vec;

    // Loop through the CIGAR string, process operations, detect SVs (primary
    // only), update clipped base support, calculate sequence identity for
    // potential duplications (primary only), and calculate
    // the clipped base support and mismatch rate
    int32_t ref_pos;
    int32_t ref_end;
    int match_count = 0;
    int mismatch_count = 0;
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
                std::string ins_seq_str = "";
                uint8_t* seq_ptr = bam_get_seq(alignment);
                for (int j = 0; j < op_len; j++) {
                    ins_seq_str += seq_nt16_str[bam_seqi(seq_ptr, query_pos + j)];
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
                std::lock_guard<std::mutex> lock(this->sv_mtx);
                sv_calls.add(chr, ref_pos, ref_end, DEL, ".", "CIGARDEL");
            }

        // Check if the CIGAR operation is a soft clip
        } else if (op == BAM_CSOFT_CLIP) {

            // Update the clipped base support
            std::lock_guard<std::mutex> lock(this->sv_mtx);
            sv_calls.updateClippedBaseSupport(chr, pos);

            // Update the query alignment start position
            if (!first_op) {
                query_start = query_pos + op_len;
                first_op = true;
            }
        }

        // Update match/mismatch counts
        if (op == BAM_CEQUAL) {
            match_count += op_len;
        } else if (op == BAM_CDIFF) {
            mismatch_count += op_len;
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

            // Compare the two sequences and update the mismatch count
            for (int j = 0; j < op_len; j++) {
                if (cmatch_seq_str[j] != cmatch_ref_str[j]) {
                    mismatch_count++;
                } else {
                    match_count++;
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

    // Calculate the mismatch rate
    double mismatch_rate = (double)mismatch_count / (double)(match_count + mismatch_count);
    std::cout << "Mismatch count: " << mismatch_count << std::endl;
    std::cout << "Match count: " << match_count << std::endl;
    std::cout << "Mismatch rate (%): " << mismatch_rate * 100 << std::endl;

    return std::tuple<double, int32_t, int32_t>(mismatch_rate, query_start, query_end);
}

// Detect SVs from split read alignments (primary and supplementary) and
// directly from the CIGAR string
SVData SVCaller::run()
{
    // Open the BAM file
    std::string bam_filepath = this->input_data->getLongReadBam();

    // Get the region data
    bool whole_genome = this->input_data->getWholeGenome();
    std::vector<std::string> regions;
    if (whole_genome) {
        regions = this->input_data->getRefGenomeChromosomes();
    } else {
        regions.push_back(this->input_data->getRegion());
    }
    int num_regions = regions.size();

    // Loop through each region and detect SVs
    std::cout << "Detecting SVs from " << num_regions << " region(s)..." << std::endl;
    int region_count = 0;
    auto start1 = std::chrono::high_resolution_clock::now();
    SVData sv_calls;
    PrimaryMap primary_alignments;
    SuppMap supplementary_alignments;
    int all_region_sv_count = 0;
    for (const auto& region : regions) {
        // std::cout << "Extracting alignments for region: " << region << std::endl;
        // Split the region into chunks and process each chunk in a separate
        // thread
        int num_threads = this->input_data->getThreadCount();

        // If a region range is specified (e.g., chr1:1-1000), use a single
        // thread
        std::vector<std::string> region_chunks;
        if (region.find(":") != std::string::npos) {
            num_threads = 1;
            region_chunks.push_back(region);
        } else {
            // Get the chromosome length
            int chr_len = this->input_data->getRefGenomeChromosomeLength(region);

            // Split the region into chunks
            int chunk_size = chr_len / num_threads;
            for (int i = 0; i < num_threads; i++) {
                int start = i * chunk_size + 1;  // 1-based
                int end = start + chunk_size;
                if (i == num_threads - 1) {
                    end = chr_len;
                }
                std::string chunk = region + ":" + std::to_string(start) + "-" + std::to_string(end);
                region_chunks.push_back(chunk);
            }
        }

        // Loop through the chunks and process each chunk in a separate thread
        std::vector<RegionData> sv_calls_vec;
        std::vector<std::thread> threads;

        // Vector of futures
        std::vector<std::future<RegionData>> futures;
        for (const auto& sub_region : region_chunks) {
            // Detect SVs from the sub-region in a separate thread using a
            // future
            std::future<RegionData> future = std::async(std::launch::async, &SVCaller::detectSVsFromRegion, this, sub_region);

            // Add the future to the list of futures
            futures.push_back(std::move(future));
        }

        // Get the SV calls from the futures
        int region_sv_count = 0;
        for (auto& future : futures) {
            // Wait for the future to be ready
            future.wait();

            // Get the SV region data from the future
            RegionData sv_calls_region = std::move(future.get());
            sv_calls_vec.push_back(sv_calls_region);

            // Update the SV count
            region_sv_count += std::get<0>(sv_calls_region).totalCalls();
            all_region_sv_count += std::get<0>(sv_calls_region).totalCalls();
        }

        // Combine the SV calls, primary alignments, and supplementary
        // alignments
        for (auto it = sv_calls_vec.begin(); it != sv_calls_vec.end(); ++it) {
            // Print the number of SVs before combining
            sv_calls.concatenate(std::get<0>(*it));

            // Combine the primary alignments
            for (auto it2 = std::get<1>(*it).begin(); it2 != std::get<1>(*it).end(); ++it2) {
                primary_alignments[it2->first] = it2->second;
            }

            // Combine the supplementary alignments
            for (auto it3 = std::get<2>(*it).begin(); it3 != std::get<2>(*it).end(); ++it3) {
                auto it4 = supplementary_alignments.find(it3->first);
                if (it4 == supplementary_alignments.end()) {
                    supplementary_alignments[it3->first] = it3->second;
                } else {
                    supplementary_alignments[it3->first].insert(supplementary_alignments[it3->first].end(), it3->second.begin(), it3->second.end());
                }
            }
        }

        // Increment the region count
        region_count++;
        std::cout << "Extracted aligments for " << region_count << " of " << num_regions << " regions." << std::endl;
    }

    // Run split-read SV detection in a single thread
    std::cout << "Detecting SVs from split-read alignments..." << std::endl;
    this->detectSVsFromSplitReads(sv_calls, primary_alignments, supplementary_alignments);

    auto end1 = std::chrono::high_resolution_clock::now();
    std::cout << "Finished detecting " << sv_calls.totalCalls() << " SVs from " << num_regions << " region(s). Elapsed time: " << getElapsedTime(start1, end1) << std::endl;

    // Return the SV calls
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

        // Get the mismatch rate
        double primary_mismatch_rate = std::get<6>(primary_alignment);

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

            // Get the mismatch rate
            double supp_mismatch_rate = std::get<6>(supp_alignment);

            // Determine if there is overlap between the primary and
            // supplementary query sequences
            int32_t overlap_start = std::max(primary_query_start, supp_query_start);
            int32_t overlap_end = std::min(primary_query_end, supp_query_end);
            int32_t overlap_length = overlap_end - overlap_start;
            if (overlap_length > 0) {
                // Choose the alignment with the lower mismatch rate,
                // and trim the overlap from the other alignment
                if (primary_mismatch_rate < supp_mismatch_rate) {
                    // Trim the overlap from the supplementary alignment
                    supp_start = supp_end - overlap_length;
                } else {
                    // Trim the overlap from the primary alignment
                    primary_end = primary_start + overlap_length;
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
