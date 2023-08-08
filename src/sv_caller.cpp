//
// sv_caller.cpp:
// Detect SVs from long read alignments
//

#include "sv_caller.h"

#include <htslib/sam.h>

#include <iostream>
#include <string>
#include <vector>
#include <map>

SVCaller::SVCaller(Common common)
{
    this->common = common;
}

// Detect SVs and return SV type by start and end position
SVMap SVCaller::run(CNVMap cnv_calls)
{
    // Test if a key exists in the CNV map
    std::pair<std::string, int> key("chr3", 60771676);
    if (cnv_calls.getCNVCalls().find(key) != cnv_calls.getCNVCalls().end()) {
        std::cout << "[2] Found key" << std::endl;
    } else {
        std::cout << "[2] Did not find key" << std::endl;
    }

    // Open the BAM file
    samFile *fp_in = sam_open(this->common.getBAMFilepath().c_str(), "r");
    if (fp_in == NULL) {
        std::cerr << "ERROR: failed to open " << this->common.getBAMFilepath() << std::endl;
        exit(1);
    }

    // Get the header
    bam_hdr_t *bamHdr = sam_hdr_read(fp_in);

    // Get the region
    std::string region = this->common.getRegion();

    // Get the index
    hts_idx_t *idx = sam_index_load(fp_in, this->common.getBAMFilepath().c_str());
    if (idx == NULL) {
        std::cerr << "ERROR: failed to load index for " << this->common.getBAMFilepath() << std::endl;
        exit(1);
    }

    // Set up the iterator
    hts_itr_t *itr = sam_itr_querys(idx, bamHdr, region.c_str());

    // Set up the read
    bam1_t *bam1 = bam_init1();

    // Set the minimum mapping quality
    int min_mapq = 20;

    // Set the minimum alignment length
    int min_alignment_length = 50;
    //int min_alignment_length = 0;

    // Iterate through the alignments in the region and get the SV calls
    SVMap sv_calls;
    int num_alignments = 0;
    int num_sv_calls = 0;
    //while (sam_read1(fp_in, bamHdr, bam1) >= 0) {
    while (sam_itr_next(fp_in, itr, bam1) >= 0) {
        
        // Skip secondary alignments
        if (bam1->core.flag & BAM_FSECONDARY) {            
            continue;
        }

        // Skip alignments with mapping quality less than min_mapq
        if (bam1->core.qual < min_mapq) {
            continue;
        }

        // Get insertion and deletion start and end positions from the CIGAR
        // string
        int32_t ins_start = -1;
        int32_t ins_end = -1;
        int32_t del_start = -1;
        int32_t del_end = -1;
        int32_t ref_pos = bam1->core.pos;
        uint32_t *cigar = bam_get_cigar(bam1);

        std::string chr = bamHdr->target_name[bam1->core.tid];
        //std::cout << "SV chr: " << chr << std::endl;
        //char *chr = bamHdr->target_name[bam1->core.tid];
        //char chr[100];
        //strcpy(chr, bamHdr->target_name[bam1->core.tid]);

        // Define the maximum distance between consecutive insertions or
        // deletions to be considered a single SV
        int max_indel_distance = 100; // 100 bp
        int ins_distance = 0;
        int prev_ins_end = -1;
        int del_distance = 0;
        int prev_del_end = -1;

        for (uint32_t i = 0; i < bam1->core.n_cigar; i++) {

            // Get the CIGAR operation
            int cigar_op = bam_cigar_op(cigar[i]);

            // Get the CIGAR operation length
            int cigar_op_len = bam_cigar_oplen(cigar[i]);

            // Check if the CIGAR operation is an insertion
            if (cigar_op == BAM_CINS) {

                // Check if the insertion start position has not been set
                if (ins_start == -1) {
                    ins_start = ref_pos;
                }

                // Set the insertion end position
                ins_end = ref_pos + cigar_op_len;

                // Get the distance between consecutive insertions
                if (prev_ins_end != -1) {
                    ins_distance = ref_pos - prev_ins_end;
                }

                // Set the previous insertion end position
                prev_ins_end = ins_end;

                // Set the insertion end position
                //ins_end = ref_pos + cigar_op_len;
                //ins_end = ins_start + cigar_op_len;
            }

            // Check if the CIGAR operation is a deletion
            else if (cigar_op == BAM_CDEL) {

                // Check if the deletion start position has not been set
                if (del_start == -1) {
                    del_start = ref_pos;
                }

                // Set the deletion end position
                del_end = ref_pos + cigar_op_len;
                //del_end = del_start + cigar_op_len;

                // Get the distance between consecutive deletions
                if (prev_del_end != -1) {
                    del_distance = ref_pos - prev_del_end;
                }

                // Set the previous deletion end position
                prev_del_end = del_end;
            }

            // Update the reference position based on the CIGAR operation and length
            if (cigar_op == BAM_CMATCH || cigar_op == BAM_CEQUAL || cigar_op == BAM_CDIFF) {
                ref_pos += cigar_op_len;
            } else if (cigar_op == BAM_CDEL || cigar_op == BAM_CREF_SKIP) {
                ref_pos += cigar_op_len;
            } else if (cigar_op == BAM_CINS || cigar_op == BAM_CSOFT_CLIP) {
                // Do nothing
            } else if (cigar_op == BAM_CHARD_CLIP) {
                // Do nothing
            } else {
                std::cerr << "ERROR: unknown CIGAR operation " << cigar_op << std::endl;
                exit(1);
            }

            // Process the insertion and deletion start and end positions

            // Check if an insertion has been found
            if (ins_start != -1 && ins_end != -1) {

                // Check if the insertion is longer than the minimum alignment length
                if (ins_end - ins_start > min_alignment_length) {

                    // Get the predicted SV type from the CNV calls
                    int sv_type = cnv_calls.getSVType(chr, ins_start, ins_end);
                    // std::cout << "cnv type: " << sv_type << std::endl;
                    // std::cout << "sv type: " << 1 << std::endl;
                    if (sv_type == 1) {
                        std::cout << "Found insertion supported by CNV calls:" << std::endl;
                        std::cout << "CIGAR INS at " << chr << ":" << ins_start << "-" << ins_end << ", length " << ins_end - ins_start << std::endl;
                    }

                    // Add the insertion to the SV calls
                    sv_calls.addSVCall(chr, ins_start, ins_end, 1);
                    num_sv_calls++;
                }

                // Reset the insertion start and end positions if the distance
                // between consecutive insertions is greater than the maximum
                // indel distance
                if (ins_distance > max_indel_distance) {
                    ins_start = -1;
                    ins_end = -1;
                    ins_distance = 0;
                } else {
                    ins_distance = ins_end - ins_start;
                }
            }

            // Check if a deletion has been found
            if (del_start != -1 && del_end != -1) {

                // Check if the deletion is longer than the minimum alignment length
                if (del_end - del_start > min_alignment_length) {

                    // Get the predicted SV type from the CNV calls
                    int sv_type = cnv_calls.getSVType(chr, del_start, del_end);
                    // std::cout << "cnv type: " << sv_type << std::endl;
                    // std::cout << "sv type: " << 2 << std::endl;
                    if (sv_type == 2) {
                        std::cout << "Found deletion supported by CNV calls" << std::endl;
                        std::cout << "CIGAR DEL at " << chr << ":" << del_start << "-" << del_end << ", length " << del_end - del_start << std::endl;
                    }

                    // Add the deletion to the SV calls
                    sv_calls.addSVCall(chr, del_start, del_end, 2);
                    num_sv_calls++;
                }

                // Reset the deletion start and end positions if the distance
                // between consecutive deletions is greater than the maximum
                // indel distance
                if (del_distance > max_indel_distance) {
                    del_start = -1;
                    del_end = -1;
                    del_distance = 0;
                }
            }
        }

        // Increment the number of alignment records processed
        num_alignments++;

        // Print the progress as a wait cursor
        if (num_alignments % 100000 == 0) {
            std::cout << "\r" << num_alignments << " alignment records processed";
            std::cout.flush();
        }
    }

    // Print the number of alignments processed
    std::cout << "\r" << num_alignments << " alignment records processed" << std::endl;

    // Destroy the iterator
    hts_itr_destroy(itr);

    // Destroy the index
    hts_idx_destroy(idx);

    // Close the BAM file
    sam_close(fp_in);

    // Destroy the alignment
    bam_destroy1(bam1);

    // Destroy the header
    bam_hdr_destroy(bamHdr);

    // Print the number of SV calls
    std::cout << num_sv_calls << " SV calls" << std::endl;

    // Return the SV calls
    return sv_calls;
}
