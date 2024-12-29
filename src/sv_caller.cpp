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
#include <condition_variable>

#include "ThreadPool.h"
#include "utils.h"
#include "sv_types.h"
#include "version.h"
#include "fasta_query.h"
/// @endcond

# define DUP_SEQSIM_THRESHOLD 0.90  // Sequence similarity threshold for duplication detection

//std::mutex bam_mutex;

int SVCaller::readNextAlignment(samFile *fp_in, hts_itr_t *itr, bam1_t *bam1)
{
    std::lock_guard<std::mutex> lock(this->shared_mutex);
    int ret = sam_itr_next(fp_in, itr, bam1);
    return ret;
}

void SVCaller::getSplitAlignments(samFile* fp_in, hts_idx_t* idx, bam_hdr_t* bamHdr, const std::string& region, std::unordered_map<std::string, GenomicRegion>& primary_map, std::unordered_map<std::string, std::vector<GenomicRegion>>& supp_map)
{
    // Create a read and iterator for the region
    bam1_t *bam1 = bam_init1();
    if (!bam1) {
        printError("ERROR: failed to initialize BAM record");
        return;
    }
    hts_itr_t *itr = sam_itr_querys(idx, bamHdr, region.c_str());
    if (!itr) {
        bam_destroy1(bam1);
        printError("ERROR: failed to query region " + region);
        return;
    }

    uint32_t primary_count = 0;
    uint32_t supplementary_count = 0;

    // Main loop to process the alignments
    uint32_t num_alignments = 0;
    while (readNextAlignment(fp_in, itr, bam1) >= 0) {

        // Skip secondary and unmapped alignments, duplicates, QC failures, and low mapping quality
        if (bam1->core.flag & BAM_FSECONDARY || bam1->core.flag & BAM_FUNMAP || bam1->core.flag & BAM_FDUP || bam1->core.flag & BAM_FQCFAIL || bam1->core.qual < this->min_mapq) {
            continue;
        }
        const std::string qname = bam_get_qname(bam1);  // Query template name

        // Process primary alignments
        if (!(bam1->core.flag & BAM_FSUPPLEMENTARY)) {
            // primary_map[qname] = itr;
            // Store chromosome (TID), start, and end positions (1-based) of the
            // primary alignment, and the strand
            primary_map[qname] = GenomicRegion{bam1->core.tid, bam1->core.pos + 1, bam_endpos(bam1), !(bam1->core.flag & BAM_FREVERSE)};
            primary_count++;

        // Process supplementary alignments
        } else if (bam1->core.flag & BAM_FSUPPLEMENTARY) {
            // supp_map[qname].push_back(itr);
            // Store chromosome (TID), start, and end positions (1-based) of the
            // supplementary alignment, and the strand
            supp_map[qname].push_back(GenomicRegion{bam1->core.tid, bam1->core.pos + 1, bam_endpos(bam1), !(bam1->core.flag & BAM_FREVERSE)});
            supplementary_count++;
        }
        num_alignments++;
    }

    // Remove primary alignments without supplementary alignments
    std::vector<std::string> to_remove;
    for (const auto& entry : primary_map) {
        const std::string& qname = entry.first;
        if (supp_map.find(qname) == supp_map.end()) {
            to_remove.push_back(qname);
        }
    }
    for (const std::string& qname : to_remove) {
        primary_map.erase(qname);
    }

    // Clean up the iterator and alignment
    hts_itr_destroy(itr);
    bam_destroy1(bam1);
    printMessage(region + ": Processed " + std::to_string(primary_map.size()) + " primary alignments with " + std::to_string(supplementary_count) + " supplementary alignments");
    // printMessage("Processed " + std::to_string(num_alignments) + " alignments with " + std::to_string(primary_count) + " primary and " + std::to_string(supplementary_count) + " supplementary alignments...");
}


void SVCaller::getAlignmentMismatchMap(samFile *fp_in, hts_idx_t *idx, bam_hdr_t *bamHdr, const GenomicRegion& region, MismatchData &mismatch_data, bool is_primary, const ReferenceGenome& ref_genome)
{
    // Create a read and iterator for the region
    bam1_t *bam1 = bam_init1();
    if (!bam1) {
        printError("ERROR: failed to initialize BAM record");
        return;
    }


    //bam_mutex.lock();
    this->shared_mutex.lock();
    hts_itr_t *itr = sam_itr_queryi(idx, region.tid, region.start - 1, region.end);
    if (!itr) {
        this->shared_mutex.unlock();
        bam_destroy1(bam1);
        printError("ERROR: failed to query region " + std::to_string(region.tid) + ":" + std::to_string(region.start) + "-" + std::to_string(region.end));
        return;
    }
    this->shared_mutex.unlock();
    //bam_mutex.unlock();


    // Find the correct alignment
    bool success = false;
    std::string fail_str = "";
    // printMessage("Looking for alignment for region: " + std::to_string(region.start) + "-" + std::to_string(region.end) + " with type: " + (is_primary ? "primary" : "supplementary") + " and strand: " + (region.strand ? "forward" : "reverse"));
    while (readNextAlignment(fp_in, itr, bam1) >= 0) {
        // Skip secondary and unmapped alignments, duplicates, QC failures, and low mapping quality
        if (bam1->core.flag & BAM_FSECONDARY || bam1->core.flag & BAM_FUNMAP || bam1->core.flag & BAM_FDUP || bam1->core.flag & BAM_FQCFAIL || bam1->core.qual < this->min_mapq) {
            continue;
        }

        // Skip if not the correct type of alignment
        if (is_primary && (bam1->core.flag & BAM_FSUPPLEMENTARY)) {
            continue;
        } else if (!is_primary && !(bam1->core.flag & BAM_FSUPPLEMENTARY)) {
            continue;
        }

        // Check the alignment start and end positions, and strand
        if (bam1->core.pos+1 == region.start && bam_endpos(bam1) == region.end && !(bam1->core.flag & BAM_FREVERSE) == region.strand) {
            // printMessage("SUCCESS: Found alignment for region: " + std::to_string(region.start) + "-" + std::to_string(region.end) + " at position: " + std::to_string(bam1->core.pos + 1) + "-" + std::to_string(bam_endpos(bam1)));
            success = true;
            break;
        } else {
            continue;
        }
    }

    // Check if the alignment was found
    if (!success) {
        printError("ERROR: Failed to find alignment for region: " + std::to_string(region.start) + "-" + std::to_string(region.end) + " with type: " + (is_primary ? "primary" : "supplementary") + " and strand: " + (region.strand ? "forward" : "reverse"));
        hts_itr_destroy(itr);
        bam_destroy1(bam1);
        return;
    }

    // Main loop to process the alignments
    std::vector<int> match_map(bam1->core.l_qseq, 0);  // Query position to match/mismatch (1/0) map
    uint32_t query_start = 0;
    uint32_t query_end = 0;
    uint32_t query_pos = 0;
    bool first_op = true;

    // Process mismatches in the CIGAR string
    const std::string chr = bamHdr->target_name[bam1->core.tid];
    hts_pos_t pos = bam1->core.pos;  // 0-based position
    uint32_t* cigar = bam_get_cigar(bam1);  // CIGAR array
    int cigar_len = bam1->core.n_cigar;
    for (int i = 0; i < cigar_len; i++) {
        int op = bam_cigar_op(cigar[i]);  // CIGAR operation
        int op_len = bam_cigar_oplen(cigar[i]);  // CIGAR operation length
        
        // Update match/mismatch query map
        int MATCH = 1;
        int MISMATCH = -1;
        if (op == BAM_CEQUAL) {
            for (int j = 0; j < op_len; j++) {
                match_map[query_pos + j] = MATCH;
            }
        } else if (op == BAM_CDIFF) {
            for (int j = 0; j < op_len; j++) {
                match_map[query_pos + j] = MISMATCH;
            }
        } else if (op == BAM_CMATCH) {
            // Get the read sequence
            uint8_t* seq_ptr = bam_get_seq(bam1);
            std::string cmatch_seq_str = "";
            for (int j = 0; j < op_len; j++) {
                cmatch_seq_str += seq_nt16_str[bam_seqi(seq_ptr, query_pos + j)];
            }

            // Get the corresponding reference sequence
            int cmatch_pos = pos + 1;  // Querying the reference genome is 1-based
            std::string cmatch_ref_str = ref_genome.query(chr, cmatch_pos, cmatch_pos + op_len - 1);

            // Check that the two sequence lengths are equal
            if (cmatch_seq_str.length() != cmatch_ref_str.length()) {
                printError("ERROR: Sequence lengths do not match for CIGAR operation: " + std::to_string(op));
                hts_itr_destroy(itr);
                bam_destroy1(bam1);
                return;
            }

            // Compare the two sequences and update the mismatch map
            for (int j = 0; j < op_len; j++) {
                if (cmatch_seq_str[j] != cmatch_ref_str[j]) {
                    try {
                        match_map.at(query_pos + j) = MISMATCH;
                    } catch (const std::out_of_range& e) {
                        printError("ERROR: Out of range exception for query position: " + std::to_string(query_pos + j) + " with read length: " + std::to_string(bam1->core.l_qseq) + " and array size: " + std::to_string(match_map.size()) + " for CIGAR operation: " + std::to_string(op) + " with length: " + std::to_string(op_len));

                        // Exit the program
                        hts_itr_destroy(itr);
                        bam_destroy1(bam1);
                        
                        return;
                    }
                    // match_map[query_pos + j] = MISMATCH;
                } else {
                    match_map[query_pos + j] = MATCH;
                }
            }
        } else if (first_op && (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP)) {
            query_start = query_pos + op_len;
            first_op = false;
        }
        
        // Update the reference position
        // https://samtools.github.io/hts-specs/SAMv1.pdf
        if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP || op == BAM_CEQUAL || op == BAM_CDIFF) {
            pos += op_len;
        }

        // Update the query position
        if (op == BAM_CMATCH || op == BAM_CINS || op == BAM_CSOFT_CLIP || op == BAM_CEQUAL || op == BAM_CDIFF) {
            query_pos += op_len;
        }
    }
    query_end = query_pos;
    
    // Clean up the iterator and alignment
    hts_itr_destroy(itr);
    bam_destroy1(bam1);

    // Update the mismatch data
    mismatch_data.query_start = query_start;
    mismatch_data.query_end = query_end;
    mismatch_data.match_map = std::move(match_map);
}


void SVCaller::detectCIGARSVs(samFile* fp_in, hts_idx_t* idx, bam_hdr_t* bamHdr, const std::string& region, std::vector<SVCall>& sv_calls, const std::vector<uint32_t>& pos_depth_map, const ReferenceGenome& ref_genome)
{
    // Create a read and iterator for the region
    bam1_t *bam1 = bam_init1();
    if (!bam1) {
        printError("ERROR: failed to initialize BAM record");
        return;
    }
    hts_itr_t *itr = sam_itr_querys(idx, bamHdr, region.c_str());
    if (!itr) {
        bam_destroy1(bam1);
        printError("ERROR: failed to query region " + region);
        return;
    }

    // Main loop to process the alignments
    while (readNextAlignment(fp_in, itr, bam1) >= 0) {

        // Skip secondary and unmapped alignments, duplicates, QC failures, and low mapping quality
        if (bam1->core.flag & BAM_FSECONDARY || bam1->core.flag & BAM_FUNMAP || bam1->core.flag & BAM_FDUP || bam1->core.flag & BAM_FQCFAIL || bam1->core.qual < this->min_mapq) {
            continue;
        }

        // Process the alignment
        bool primary = !(bam1->core.flag & BAM_FSUPPLEMENTARY);
        this->detectSVsFromCIGAR(bamHdr, bam1, sv_calls, primary, pos_depth_map, ref_genome);
    }

    // Clean up the iterator and alignment
    hts_itr_destroy(itr);
    bam_destroy1(bam1);
}

double SVCaller::calculateMismatchRate(const MismatchData& mismatch_data)
{
    int start = mismatch_data.query_start;
    int end = mismatch_data.query_end;
    const std::vector<int>& mismatch_map = mismatch_data.match_map;
    start = std::max(start, 0);
    end = std::min(end, (int32_t)mismatch_map.size() - 1);
    int match_count = 0;
    int mismatch_count = 0;
    int MATCH = 1;
    int MISMATCH = -1;
    for (int i = start; i <= end; i++) {
        if (mismatch_map[i] == MATCH) {
            match_count++;
        } else if (mismatch_map[i] == MISMATCH) {
            mismatch_count++;
        }
    }

    // Avoid division by zero
    if (match_count + mismatch_count == 0) {
        return 0.0;
    }

    double mismatch_rate = static_cast<double>(mismatch_count) / static_cast<double>(match_count + mismatch_count);

    return mismatch_rate;
}

void SVCaller::detectSVsFromCIGAR(bam_hdr_t* header, bam1_t* alignment, std::vector<SVCall>& sv_calls, bool is_primary, const std::vector<uint32_t>& pos_depth_map, const ReferenceGenome& ref_genome)
{
    std::string chr = header->target_name[alignment->core.tid];  // Chromosome name
    uint32_t pos = (uint32_t)alignment->core.pos;  // Leftmost position of the alignment in the reference genome (0-based)
    uint32_t* cigar = bam_get_cigar(alignment);  // CIGAR array
    int cigar_len = alignment->core.n_cigar;
    uint32_t query_pos = 0;

    // Loop through the CIGAR string, process operations, detect SVs (primary
    // only), and calculate sequence identity for potential duplications (primary only)
    uint32_t ref_pos;
    uint32_t ref_end;
    double default_lh = 0.0;
    // List of ambiguous bases
    const std::string amb_bases = "RYKMSWBDHV";
    for (int i = 0; i < cigar_len; i++) {

        int op = bam_cigar_op(cigar[i]);  // CIGAR operation
        int op_len = bam_cigar_oplen(cigar[i]);  // CIGAR operation length

        if (op_len == 0) {
            printError("Warning: Encountered CIGAR operation with length 0 at position " + std::to_string(pos+1) + " in chromosome " + chr);
            continue;
        }
        
        // Process the CIGAR operation
        if (op == BAM_CINS && is_primary) {

            // Get the sequence of the insertion from the query
            std::string ins_seq_str(op_len, ' ');
            for (int j = 0; j < op_len; j++) {
                // Replace ambiguous bases with N
                char base = seq_nt16_str[bam_seqi(bam_get_seq(alignment), query_pos + j)];
                if (amb_bases.find(base) != std::string::npos) {
                    ins_seq_str[j] = 'N';
                } else {
                    ins_seq_str[j] = base;
                }
                // Get the sequence character from the query
                // ins_seq_str[j] = seq_nt16_str[bam_seqi(bam_get_seq(alignment), query_pos + j)];
            }

            // To determine whether the insertion is a duplication, check
            // for sequence identity between the insertion and the
            // reference genome (duplications are typically >= 90%):
            // Loop through the reference sequence and calculate the
            // sequence identity +/- insertion length from the insertion
            // position.
            // bool is_duplication = false;
            // int ins_ref_pos;
            // uint32_t dup_start = std::max(0, (int)pos - op_len);
            // for (uint32_t j = dup_start; j <= pos; j++) {

            //     // Get the string for the window (1-based coordinates)
            //     ins_ref_pos = j + 1;
            //     std::string window_str = ref_genome.query(chr, ins_ref_pos, ins_ref_pos + op_len - 1);

            //     // Continue if the window string is empty (out-of-range)
            //     if (window_str == "") {
            //         continue;
            //     }

            //     // Calculate the sequence identity
            //     int num_matches = 0;
            //     for (int k = 0; k < op_len; k++) {
            //         if (ins_seq_str[k] == window_str[k] && ins_seq_str[k] != 'N' && window_str[k] != 'N') {
            //             num_matches++;
            //         }
            //     }
            //     float seq_identity = (float)num_matches / (float)op_len;

            //     // Check if the target sequence identity is reached
            //     if (seq_identity >= DUP_SEQSIM_THRESHOLD) {
            //         is_duplication = true;
            //         break;
            //     }
            // }

            // Calculate the sequence identity at the insertion position +/-
            // length
            // Before the insertion
            if (pos >= (uint32_t)op_len-1)
            {
                uint32_t bp1 = pos - (op_len - 1);
                uint32_t bp2 = pos;
                const std::string& window_str = ref_genome.query(chr, bp1 + 1, bp2 + 1);
                if (window_str.length() > 0)
                {
                    int num_matches = 0;
                    for (int k = 0; k < op_len; k++)
                    {
                        if (ins_seq_str[k] == window_str[k] && ins_seq_str[k] != 'N' && window_str[k] != 'N')
                        {
                            num_matches++;
                        }
                    }
                    float seq_identity = (float)num_matches / (float)op_len;
                    if (seq_identity >= DUP_SEQSIM_THRESHOLD)
                    {
                        uint32_t dup_bp1 = bp1 + 1;
                        uint32_t dup_bp2 = std::min(dup_bp1 + op_len - 1, ref_genome.getChromosomeLength(chr));
                        int read_depth = this->calculateReadDepth(pos_depth_map, dup_bp1, dup_bp2);
                        addSVCall(sv_calls, dup_bp1, dup_bp2, "DUP", ins_seq_str, "LSEQSIM", "./.", default_lh, read_depth);

                        // Continue to the next CIGAR operation
                        continue;
                    }
                }
            }

            // After the insertion
            if (pos + op_len < ref_genome.getChromosomeLength(chr))
            {
                uint32_t bp1 = pos + 1;
                uint32_t bp2 = bp1 + op_len - 1;
                const std::string& window_str = ref_genome.query(chr, bp1 + 1, bp2 + 1);
                if (window_str.length() > 0)
                {
                    int num_matches = 0;
                    for (int k = 0; k < op_len; k++)
                    {
                        if (ins_seq_str[k] == window_str[k] && ins_seq_str[k] != 'N' && window_str[k] != 'N')
                        {
                            num_matches++;
                        }
                    }
                    float seq_identity = (float)num_matches / (float)op_len;
                    if (seq_identity >= DUP_SEQSIM_THRESHOLD)
                    {
                        uint32_t dup_bp1 = bp1 + 1;
                        uint32_t dup_bp2 = std::min(dup_bp1 + op_len - 1, ref_genome.getChromosomeLength(chr));
                        int read_depth = this->calculateReadDepth(pos_depth_map, dup_bp1, dup_bp2);
                        addSVCall(sv_calls, dup_bp1, dup_bp2, "DUP", ins_seq_str, "RSEQSIM", "./.", default_lh, read_depth);

                        // Continue to the next CIGAR operation
                        continue;
                    }
                }
            }

            // Add as an insertion
            // For read depth calculation, use the previous and current
            // positions (1-based)
            int read_depth = this->calculateReadDepth(pos_depth_map, std::max(1, (int)pos), pos + 1);
            uint32_t ins_pos = pos + 1;
            uint32_t ins_end = ins_pos + op_len - 1;
            addSVCall(sv_calls, ins_pos, ins_end, "INS", ins_seq_str, "CIGARINS", "./.", default_lh, read_depth);

            // Determine whether to use a symbolic allele (>50bp) or the
            // actual sequence
            // if (op_len > 50) {
            //     ins_seq_str = "<INS>";
            // } else {
            //     ins_seq_str = ins_seq_str;
            // }

            // Add to SV calls (1-based) with the appropriate SV type
            // ref_pos = pos+1;

            // // For insertions, the reference end position is the same as the
            // // reference position
            // // For duplications, the reference end position is the same as
            // // the reference position plus the length of the insertion
            // ref_end = ref_pos + op_len - 1;
            // if (is_duplication) {
            //     uint32_t bp1 = ref_pos;
            //     uint32_t bp2 = std::min(ref_pos + op_len - 1, ref_genome.getChromosomeLength(chr));
            //     int read_depth = this->calculateReadDepth(pos_depth_map, bp1, bp2);
            //     addSVCall(sv_calls, ref_pos, ref_end, "DUP", ins_seq_str, "CIGARDUP", "./.", default_lh, read_depth);
            // } else {
            //     uint32_t bp1 = std::max(1, (int)ref_pos - 1);
            //     uint32_t bp2 = ref_pos;
            //     int read_depth = this->calculateReadDepth(pos_depth_map, bp1, bp2);
            //     addSVCall(sv_calls, ref_pos, ref_end, "INS", ins_seq_str, "CIGARINS", "./.", default_lh, read_depth);
            // }

        // Check if the CIGAR operation is a deletion
        } else if (op == BAM_CDEL && is_primary) {

            ref_pos = pos+1;
            ref_end = ref_pos + op_len -1;
            // printMessage("Test2");
            int read_depth = this->calculateReadDepth(pos_depth_map, ref_pos, ref_end);
            addSVCall(sv_calls, ref_pos, ref_end, "DEL", "<DEL>", "CIGARDEL", "./.", default_lh, read_depth);
        }

        // Update the reference position
        // https://samtools.github.io/hts-specs/SAMv1.pdf
        if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP || op == BAM_CEQUAL || op == BAM_CDIFF) {
            pos += op_len;
        }
        
        // Update the query position
        if (op == BAM_CMATCH || op == BAM_CINS || op == BAM_CSOFT_CLIP || op == BAM_CEQUAL || op == BAM_CDIFF) {
            query_pos += op_len;
        }
    }
}

void SVCaller::processChromosome(const std::string& chr, const CHMM& hmm, std::vector<SVCall>& chr_sv_calls, const InputData& input_data, const ReferenceGenome& ref_genome)
{
    // int filter_threshold = 4;  // Minimum number of supporting reads for an
    // SV call
    int filter_threshold = 10;  // Minimum number of supporting reads for an SV call
    bool single_chr = input_data.getChromosome() != "";

    // Open the BAM file
    std::string bam_filepath = input_data.getLongReadBam();
    samFile *fp_in = sam_open(bam_filepath.c_str(), "r");
    if (!fp_in) {
        printError("ERROR: failed to open " + bam_filepath);
        return;
    }

    // Load the header
    bam_hdr_t *bamHdr = sam_hdr_read(fp_in);
    if (!bamHdr) {
        sam_close(fp_in);
        printError("ERROR: failed to read header from " + bam_filepath);
        return;
    }

    // Load the index
    hts_idx_t *idx = sam_index_load(fp_in, bam_filepath.c_str());
    if (!idx) {
        bam_hdr_destroy(bamHdr);
        sam_close(fp_in);
        printError("ERROR: failed to load index for " + bam_filepath);
        return;
    }
    BamFileGuard bam_guard(fp_in, idx, bamHdr);  // Guard to close the BAM file

    // Set the region to process
    std::string region = chr;
    uint32_t chr_len = ref_genome.getChromosomeLength(chr);
    if (input_data.isRegionSet()) {

        // Use one chunk for the specified region
        std::pair<int32_t, int32_t> region_data = input_data.getRegion();
        int region_start = region_data.first;
        int region_end = region_data.second;
        region = chr + ":" + std::to_string(region_start) + "-" + std::to_string(region_end);
        
    }

    // Load chromosome data for copy number predictions
    printMessage(chr + ": Loading chromosome data...");
    CNVCaller cnv_caller(this->shared_mutex);
    std::vector<uint32_t> chr_pos_depth_map(chr_len+1, 0);  // 1-based index
    int thread_count = input_data.getThreadCount();
    double mean_chr_cov = cnv_caller.calculateMeanChromosomeCoverage(chr, chr_pos_depth_map, bam_filepath, thread_count, single_chr);
    if (mean_chr_cov == 0.0 || chr_pos_depth_map.size() == 0) {
        return;
    }

    // Detect SVs from the CIGAR strings
    printMessage(chr + ": CIGAR SVs...");
    this->detectCIGARSVs(fp_in, idx, bamHdr, region, chr_sv_calls, chr_pos_depth_map, ref_genome);

    printMessage(chr + ": Merging CIGAR...");
    filterSVsWithLowSupport(chr_sv_calls, filter_threshold);
    mergeSVs(chr_sv_calls);
    int region_sv_count = getSVCount(chr_sv_calls);
    printMessage("Total SVs detected from CIGAR string: " + std::to_string(region_sv_count));

    // Testing on HG002 whole genome
    // Run copy number variant predictions on the SVs detected from the
    // CIGAR string, using a minimum CNV length threshold
    if (region_sv_count > 0) {
        printMessage(chr + ": CIGAR predictions...");
        cnv_caller.runCIGARCopyNumberPrediction(chr, chr_sv_calls, hmm, mean_chr_cov, chr_pos_depth_map, input_data);
    }

    // Run split-read SV and copy number variant predictions
    printMessage(chr + ": Split read SVs...");
    this->detectSVsFromSplitReads(region, fp_in, idx, bamHdr, chr_sv_calls, cnv_caller, hmm, mean_chr_cov, chr_pos_depth_map, input_data, ref_genome);

    // Merge the SV calls from the current region
    printMessage(chr + ": Merging split reads...");
    filterSVsWithLowSupport(chr_sv_calls, filter_threshold);
    mergeSVs(chr_sv_calls);

    // Run a final merge on the combined SV calls
    printMessage(chr + ": Merging final calls...");
    mergeSVs(chr_sv_calls);
    printMessage("Completed chromosome " + chr);
}

void SVCaller::run(const InputData& input_data)
{
    // Set up the reference genome
    printMessage("Loading the reference genome...");
    const std::string ref_filepath = input_data.getRefGenome();
    ReferenceGenome ref_genome(this->shared_mutex);
    ref_genome.setFilepath(ref_filepath);

    // Get the chromosomes
    std::vector<std::string> chromosomes;
    if (input_data.isSingleChr()) {
        chromosomes.push_back(input_data.getChromosome());
    } else {
        chromosomes = ref_genome.getChromosomes();
    }
    
    // Read the HMM from the file
    std::string hmm_filepath = input_data.getHMMFilepath();
    std::cout << "Reading HMM from file: " << hmm_filepath << std::endl;
    const CHMM& hmm = ReadCHMM(hmm_filepath.c_str());

    // Use multi-threading across chromosomes unless a single chromosome is
    // specified
    int max_threads = 1;
    if (!input_data.isSingleChr()) {
        max_threads = input_data.getThreadCount();
        std::cout << "Using " << max_threads << " threads for processing..." << std::endl;
    }
    ThreadPool pool(max_threads);

    // Shared resources
    std::unordered_map<std::string, std::vector<SVCall>> whole_genome_sv_calls;
    //std::mutex sv_mutex;
    //std::mutex snp_mutex;
    //std::mutex pfb_mutex;

    // Lambda to process a chromosome
    auto process_chr = [&](const std::string& chr) {
        try {
            std::vector<SVCall> sv_calls;
            InputData chr_input_data = input_data;  // Use a thread-local copy
            this->processChromosome(chr, hmm, sv_calls, chr_input_data, ref_genome);
            {
                //std::lock_guard<std::mutex> lock(sv_mutex);
                std::lock_guard<std::mutex> lock(this->shared_mutex);
                whole_genome_sv_calls[chr] = std::move(sv_calls);
            }
            // printMessage("Completed chromosome " + chr);
        } catch (const std::exception& e) {
            printError("Error processing chromosome " + chr + ": " + e.what());
        } catch (...) {
            printError("Unknown error processing chromosome " + chr);
        }
    };

    // Submit tasks to the thread pool and track futures
    std::vector<std::future<void>> futures;
    for (const auto& chr : chromosomes) {
        futures.emplace_back(pool.enqueue([&, chr] {
            printMessage("Processing chromosome " + chr);
            process_chr(chr);
        }));
    }

    // Wait for all tasks to complete
    for (auto& future : futures) {
        try {
            future.get();
            printMessage("Chromosome task completed.");
        } catch (const std::exception& e) {
            printError("Error processing chromosome task: " + std::string(e.what()));
        } catch (...) {
            printError("Unknown error processing chromosome task.");
        }
    }
    printMessage("All tasks have finished.");

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
    const std::string output_dir = input_data.getOutputDir();
    this->saveToVCF(whole_genome_sv_calls, output_dir, ref_genome);
}


// Detect SVs from split read alignments
void SVCaller::detectSVsFromSplitReads(const std::string& region, samFile* fp_in, hts_idx_t* idx, bam_hdr_t* bamHdr, std::vector<SVCall>& sv_calls, const CNVCaller& cnv_caller, const CHMM& hmm, double mean_chr_cov, const std::vector<uint32_t>& pos_depth_map, const InputData& input_data, const ReferenceGenome& ref_genome)
{
    printMessage(region + ": Getting split alignments...");
    std::unordered_map<std::string, GenomicRegion> primary_map;
    std::unordered_map<std::string, std::vector<GenomicRegion>> supp_map;
    this->getSplitAlignments(fp_in, idx, bamHdr, region, primary_map, supp_map);

    // Find split-read SV evidence
    printMessage(region + ": Finding split-read SVs...");
    int sv_count = 0;
    int current_primary = 0;
    //int primary_count = primary_map.size();
    uint32_t min_cnv_length = input_data.getMinCNVLength();
    for (auto& entry : primary_map) {
        current_primary++;
        const std::string& qname = entry.first;
        GenomicRegion& primary_region = entry.second;

        // Skip primary alignments that do not have supplementary alignments
        // if (supp_map.find(qname) == supp_map.end()) {
        //     continue;
        // }

        // Get the read match/mismatch map
        // printMessage(region + ": Getting mismatch map for " + std::to_string(current_primary) + " of " + std::to_string(primary_count) + " primary alignments...");
        //MismatchData primary_mismatches;
        //this->getAlignmentMismatchMap(fp_in, idx, bamHdr, primary_region, primary_mismatches, true, ref_genome);
        
      	// Find the largest supplementary alignment
        GenomicRegion largest_supp_region = supp_map[qname][0];
        uint32_t largest_supp_length = 0;

        // printMessage(region + ": Processing supplementary alignments for " + std::to_string(current_primary) + " of " + std::to_string(primary_count) + " primary alignments...");
        const std::string& primary_chr = bamHdr->target_name[primary_region.tid];
        for (auto it = supp_map[qname].begin(); it != supp_map[qname].end(); ++it) {
            GenomicRegion& supp_region = *it;

            // Skip if not on the primary chromosome
            if (primary_region.tid != supp_region.tid) {
                continue;
            }

            // Get the supplementary alignment information
            uint32_t supp_start = (uint32_t) supp_region.start;
            uint32_t supp_end = (uint32_t) supp_region.end;
            uint32_t supp_length = supp_end - supp_start + 1;
            if (supp_length > largest_supp_length) {
                largest_supp_length = supp_length;
                largest_supp_region = *it;
            }

            // Inversion detection
            bool is_opposite_strand = primary_region.strand != supp_region.strand;
            if (is_opposite_strand) {
                if (supp_length >= min_cnv_length) {

                    // Print error if the start position is greater than the end
                    // position
                    if (supp_start > supp_end) {
                        printError("ERROR: Invalid inversion coordinates: " + primary_chr + ":" + std::to_string(supp_start) + "-" + std::to_string(supp_end));
                        continue;
                    }

                    // printMessage(region + ": Running copy number prediction for inversion (position: " + std::to_string(supp_start) + "-" + std::to_string(supp_end) + ")...");
                    std::tuple<double, SVType, std::string, bool> result = cnv_caller.runCopyNumberPrediction(primary_chr, hmm, supp_start, supp_end, mean_chr_cov, pos_depth_map, input_data);
                    if (std::get<1>(result) == SVType::UNKNOWN) {
                        continue;
                    }

                    double supp_lh = std::get<0>(result);
                    SVType supp_type = std::get<1>(result);
                    // printMessage("Test3");
                    int read_depth = this->calculateReadDepth(pos_depth_map, supp_start, supp_end);
                    if (supp_type == SVType::NEUTRAL) {
                        addSVCall(sv_calls, supp_start, supp_end, "INV", "<INV>", "HMM", "./.", supp_lh, read_depth);
                        
                        sv_count++;
                    } else if (supp_type == SVType::DUP) {
                        addSVCall(sv_calls, supp_start, supp_end, "INVDUP", "<INV>", "HMM", "./.", supp_lh, read_depth);
                    }
                }
                // } else {
                //     // Add the inversion without running copy number predictions
                //     // (too small for predictions)
                //     // printMessage("Test4");
                //     int read_depth = this->calculateReadDepth(pos_depth_map, supp_start, supp_end);
                //     addSVCall(sv_calls, supp_start, supp_end, "INV", "<INV>", "REV", "./.", 0.0, read_depth);
                // }
            }
        }

        // Trim overlapping alignments
        //MismatchData supp_mismatches;
        //printMessage(region + ": Getting mismatch map for supplementary alignments...");
        //this->getAlignmentMismatchMap(fp_in, idx, bamHdr, largest_supp_region, supp_mismatches, false, ref_genome);

        // printMessage(region + ": Trimming overlapping alignments...");
        //trimOverlappingAlignments(primary_region, largest_supp_region, primary_mismatches, supp_mismatches);
        bool gap_exists = false;
        uint32_t boundary_left, boundary_right, gap_left, gap_right;
        if (primary_region.start < largest_supp_region.start) {  // Primary before supp
            boundary_left = primary_region.start;
            boundary_right = std::max(primary_region.end, largest_supp_region.end);
            gap_left = primary_region.end;
            gap_right = largest_supp_region.start;
            gap_exists = gap_left < gap_right;
        } else {
            boundary_left = largest_supp_region.start;
            boundary_right = std::max(primary_region.end, largest_supp_region.end);
            gap_left = largest_supp_region.end;
            gap_right = primary_region.start;
            gap_exists = gap_left < gap_right;
        }
        
        // Run copy number variant predictions on the boundary if large enough
        if (boundary_right - boundary_left >= min_cnv_length) {

            // Print error if the start position is greater than the end
            // position
            if (boundary_left > boundary_right) {
                printError("ERROR: Invalid boundary coordinates: " + primary_chr + ":" + std::to_string(boundary_left) + "-" + std::to_string(boundary_right));
                continue;
            }

            // printMessage(region + ": Running copy number prediction for boundary...");
            std::tuple<double, SVType, std::string, bool> bd_result = cnv_caller.runCopyNumberPrediction(primary_chr, hmm, boundary_left, boundary_right, mean_chr_cov, pos_depth_map, input_data);
            if (std::get<1>(bd_result) == SVType::UNKNOWN) {
                continue;
            }
            double bd_lh = std::get<0>(bd_result);
            SVType bd_type = std::get<1>(bd_result);

            // Run copy number variant predictions on the gap if it exists
            if (gap_exists && gap_right - gap_left >= min_cnv_length) {

                // Print error if the start position is greater than the end
                // position
                if (gap_left > gap_right) {
                    printError("ERROR: Invalid gap coordinates: " + primary_chr + ":" + std::to_string(gap_left) + "-" + std::to_string(gap_right));
                    continue;
                }

                // printMessage(region + ": Running copy number prediction for gap...");
                std::tuple<double, SVType, std::string, bool> gap_result = cnv_caller.runCopyNumberPrediction(primary_chr, hmm, gap_left, gap_right, mean_chr_cov, pos_depth_map, input_data);
                if (std::get<1>(gap_result) == SVType::UNKNOWN) {
                    continue;
                }
                double gap_lh = std::get<0>(gap_result);
                SVType gap_type = std::get<1>(gap_result);

                // If higher likelihood than the boundary, add the gap as the SV call
                if (gap_lh > bd_lh) {
                    // printMessage("Test5");
                    int read_depth = this->calculateReadDepth(pos_depth_map, gap_left, gap_right);
                    std::string alt_allele = gap_type == SVType::NEUTRAL ? "." : "<" + getSVTypeString(gap_type) + ">";
                    addSVCall(sv_calls, gap_left, gap_right, getSVTypeString(gap_type), alt_allele, "GAP", "./.", gap_lh, read_depth);
                } else {
                    // Add the boundary as the SV call
                    // printMessage("Test6");
                    int read_depth = this->calculateReadDepth(pos_depth_map, boundary_left, boundary_right);
                    std::string alt_allele = bd_type == SVType::NEUTRAL ? "." : "<" + getSVTypeString(bd_type) + ">";
                    addSVCall(sv_calls, boundary_left, boundary_right, getSVTypeString(bd_type), alt_allele, "BOUNDARY", "./.", bd_lh, read_depth);
                }
            } else {
                // Add the boundary as the SV call
                // printMessage("Test7");
                int read_depth = this->calculateReadDepth(pos_depth_map, boundary_left, boundary_right);
                std::string alt_allele = bd_type == SVType::NEUTRAL ? "." : "<" + getSVTypeString(bd_type) + ">";
                addSVCall(sv_calls, boundary_left, boundary_right, getSVTypeString(bd_type), alt_allele, "BOUNDARY", "./.", bd_lh, read_depth);
            }
        }
    }
}

void SVCaller::saveToVCF(const std::unordered_map<std::string, std::vector<SVCall>>& sv_calls, const std::string& output_dir, const ReferenceGenome& ref_genome) const
{
    std::cout << "Creating VCF writer..." << std::endl;
    std::string output_vcf = output_dir + "/output.vcf";
    std::cout << "Writing VCF file to " << output_vcf << std::endl;
	std::ofstream vcf_stream(output_vcf);
    if (!vcf_stream.is_open()) {
        printError("Failed to open VCF file for writing.");
        return;
    }
    
    std::string sample_name = "SAMPLE";

    std::cout << "Getting reference genome filepath..." << std::endl;
    try {
        std::string ref_fp = ref_genome.getFilepath();
        std::cout << "Reference genome filepath: " << ref_fp << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return;
    }

    // Set the header lines
    std::cout << "Getting reference genome header..." << std::endl;
    const std::string contig_header = ref_genome.getContigHeader();
    std::vector<std::string> header_lines = {
        std::string("##reference=") + ref_genome.getFilepath(),
        contig_header,
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">",
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">",
        "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">",
        "##INFO=<ID=SVMETHOD,Number=1,Type=String,Description=\"Method used to call the structural variant\">",
        "##INFO=<ID=ALN,Number=1,Type=String,Description=\"Feature used to identify the structural variant\">",
        "##INFO=<ID=HMM,Number=1,Type=Float,Description=\"HMM likelihood\">",
        "##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description=\"Number of reads supporting the variant\">",
        "##FILTER=<ID=PASS,Description=\"All filters passed\">",
        "##FILTER=<ID=LowQual,Description=\"Low quality\">",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth at the variant site (sum of start and end positions)\">",
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
    std::string sv_method = "ContextSV" + std::string(VERSION);
    std::string source = "##source=" + sv_method;
    vcf_stream << source << std::endl;

    // Loop over the header metadata lines
    for (const auto &line : header_lines) {
        vcf_stream << line << std::endl;
    }

    // Add the header line
    std::string header_line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE";
    vcf_stream << header_line << std::endl;
    std::cout << "Saving SV calls to " << output_vcf << std::endl;
    int skip_count = 0;
    int total_count = 0;
    for (const auto& pair : sv_calls) {
        std::string chr = pair.first;
        const std::vector<SVCall>& sv_calls = pair.second;
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
            int read_depth = sv_call.read_depth;
            std::string ref_allele = ".";
            int support = sv_call.support;

            // If the SV type is unknown, skip it
            if (sv_type_str == "UNKNOWN" || sv_type_str == "NEUTRAL") {
                skip_count += 1;
                continue;
            } else {
                total_count += 1;
            }

            // Deletion
            if (sv_type_str == "DEL") {
                // Get the deleted sequence from the reference genome, also including the preceding base
                int64_t preceding_pos = (int64_t) std::max(1, (int) start-1);  // Make sure the position is not negative
                // ref_allele = ref_genome.query(chr, preceding_pos, end);
                // ref_allele = this->input_data.queryRefGenome(chr,
                // preceding_pos, end);
                ref_allele = ref_genome.query(chr, preceding_pos, end);

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
                // ref_allele = this->input_data.queryRefGenome(chr,
                // preceding_pos, preceding_pos);
                ref_allele = ref_genome.query(chr, preceding_pos, preceding_pos);

                // Update the start position to the preceding base
                start = preceding_pos;

                // Update the end position to the same base for duplications and insertions
                if (sv_type_str == "DUP" || sv_type_str == "INS") {
                    end = start;
                }

                if (sv_type_str == "INS") {
                    // Check if in symbolic form
                    if (alt_allele != "<INS>") {
                        // Use the insertion sequence as the alternate allele
                        alt_allele.insert(0, ref_allele);
                    }
                    // start = preceding_pos;  // Update the position to the preceding base

                    // // Update the end position to the start position to change from
                    // // query to reference coordinates for insertions
                    // end = start;
                }
            }

            // Create the VCF parameter strings
            std::string info_str = "END=" + std::to_string(end) + ";SVTYPE=" + sv_type_str + \
                ";SVLEN=" + std::to_string(sv_length) + ";SVMETHOD=" + sv_method + ";ALN=" + data_type_str + \
                ";HMM=" + std::to_string(hmm_likelihood) + ";SUPPORT=" + std::to_string(support);
                
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

void SVCaller::trimOverlappingAlignments(GenomicRegion& primary_alignment, GenomicRegion& supp_alignment, const MismatchData& primary_mismatches, const MismatchData& supp_mismatches)
{

    // Check for overlapping read alignments
    if (primary_mismatches.query_start < supp_mismatches.query_start) {
        // Primary before supplementary in the query

        // if (primary_query_end >= supp_query_start) {
        if (primary_mismatches.query_end >= supp_mismatches.query_start) {
            // Calculate the mismatch rates at the overlapping region
            double primary_mismatch_rate = this->calculateMismatchRate(primary_mismatches);
            double supp_mismatch_rate = this->calculateMismatchRate(supp_mismatches);
            hts_pos_t overlap_length = primary_mismatches.query_end - supp_mismatches.query_start + 1;

            // Trim the ailgnment with the higher mismatch rate
            if (primary_mismatch_rate > supp_mismatch_rate) {
                // Trim the end of the primary alignment, ensuring that the new
                // end is not less than the start
                if (primary_alignment.end > overlap_length && (primary_alignment.end - overlap_length) > primary_alignment.start) {
                    // Trim the end of the primary alignment
                    primary_alignment.end = primary_alignment.end - overlap_length;
                }
            } else {
                // Trim the beginning of the supplementary alignment, ensuring
                // that the new start is not greater than the end
                if (supp_alignment.start + overlap_length < supp_alignment.end) {
                    // Trim the beginning of the supplementary alignment
                    supp_alignment.start = supp_alignment.start + overlap_length;
                }
            }
        }

    } else {
        // Supplementary before primary in the query
        if (primary_mismatches.query_start <= supp_mismatches.query_end) {
            // Calculate the mismatch rates at the overlapping region
            double primary_mismatch_rate = this->calculateMismatchRate(primary_mismatches);
            double supp_mismatch_rate = this->calculateMismatchRate(supp_mismatches);
            hts_pos_t overlap_length = supp_mismatches.query_end - primary_mismatches.query_start + 1;

            // Trim the ailgnment with the higher mismatch rate
            if (supp_mismatch_rate > primary_mismatch_rate) {
                // Trim the end of the supplementary alignment, ensuring that
                // the new end is not less than the start
                if (supp_alignment.end > overlap_length && (supp_alignment.end - overlap_length) > supp_alignment.start) {
                    // Trim the end of the supplementary alignment
                    supp_alignment.end = supp_alignment.end - overlap_length;
                }
            } else {
                // Trim the beginning of the primary alignment, ensuring that
                // the new start is not greater than the end
                if (primary_alignment.start + overlap_length < primary_alignment.end) {
                    // Trim the beginning of the primary alignment
                    primary_alignment.start = primary_alignment.start + overlap_length;
                }
            }
        }
    }
}

int SVCaller::calculateReadDepth(const std::vector<uint32_t>& pos_depth_map, uint32_t start, uint32_t end)
{
    int read_depth = 0;
    try {
        // printMessage("Read depth at start: " + std::to_string(pos_depth_map.at(start)) + " for SV at " + std::to_string(start) + "-" + std::to_string(end) + " with length " + std::to_string(end-start));
        read_depth += pos_depth_map.at(start);
    } catch (const std::out_of_range& e) {
        // std::cerr << "Warning: Start position " << start << " not found in
        // depth map." << std::endl;
        printError("Error: Start position " + std::to_string(start) + " not found in depth map.");
    }
    try {
        // printMessage("Read depth at end: " + std::to_string(pos_depth_map.at(end)) + " for SV at " + std::to_string(start) + "-" + std::to_string(end) + " with length " + std::to_string(end-start));
        read_depth += pos_depth_map.at(end);
    } catch (const std::out_of_range& e) {
        printError("Error: End position " + std::to_string(end) + " not found in depth map.");
        // std::cerr << "Warning: End position " << end << " not found in depth map of size " << pos_depth_map.size() << "." << std::endl;
    }
    // printMessage("Read depth for SV at " + std::to_string(start) + "-" + std::to_string(end) + " with length " + std::to_string(end-start) + ": " + std::to_string(read_depth));
    return read_depth;
}
