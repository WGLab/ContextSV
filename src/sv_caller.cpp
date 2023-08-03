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
SVMap SVCaller::run()
{
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
        char *chr = bamHdr->target_name[bam1->core.tid];

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
            }

            // Check if the CIGAR operation is a deletion
            if (cigar_op == BAM_CDEL) {

                // Check if the deletion start position has not been set
                if (del_start == -1) {
                    del_start = ref_pos;
                }

                // Set the deletion end position
                del_end = ref_pos + cigar_op_len;
            }

            // Check if the CIGAR operation is a match or mismatch
            if (cigar_op == BAM_CMATCH || cigar_op == BAM_CEQUAL || cigar_op == BAM_CDIFF) {

                // Increment the reference position
                ref_pos += cigar_op_len;
            }

            // Check if the CIGAR operation is a soft clip
            if (cigar_op == BAM_CSOFT_CLIP) {

                // Increment the reference position
                ref_pos += cigar_op_len;
            }

            // Check if an insertion has been found
            if (ins_start != -1 && ins_end != -1) {

                // Check if the insertion is longer than the minimum alignment length
                if (ins_end - ins_start > min_alignment_length) {

                    // Add the insertion to the SV calls
                    sv_calls.addSVCall(chr, ins_start, ins_end, 1);
                    num_sv_calls++;
                }

                // Reset the insertion start and end positions
                ins_start = -1;
                ins_end = -1;
            }

            // Check if a deletion has been found
            if (del_start != -1 && del_end != -1) {

                // Check if the deletion is longer than the minimum alignment length
                if (del_end - del_start > min_alignment_length) {

                    // Add the deletion to the SV calls
                    sv_calls.addSVCall(chr, del_start, del_end, 2);
                    num_sv_calls++;
                }

                // Reset the deletion start and end positions
                del_start = -1;
                del_end = -1;
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
