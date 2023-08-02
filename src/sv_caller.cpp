//
// sv_caller.cpp:
// Detect SVs from long read alignments
//

#include "sv_caller.h"
#include "common.h"

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
std::map<std::pair<char *, int>, std::pair<int, int>> SVCaller::run()
{
    // Open the BAM file
    samFile *fp_in = sam_open(this->common.getBAMFilepath().c_str(), "r");
    if (fp_in == NULL) {
        std::cerr << "ERROR: failed to open " << this->common.getBAMFilepath() << std::endl;
        exit(1);
    }

    // Get the header
    bam_hdr_t *bamHdr = sam_hdr_read(fp_in);

    // Set up the read
    bam1_t *bam1 = bam_init1();

    // Variables to check if the data is sorted
    int32_t last_tid = -1;
    int32_t last_pos = -1;
    bool sorted = true;

    // SV calls by start and end position (key=[chromosome, SV start position],
    // value=[SV end position, SV type])
    std::map<std::pair<char *, int>, std::pair<int, int>> sv_calls;

    // Set the minimum mapping quality
    int min_mapq = 20;

    // Set the minimum alignment length
    int min_alignment_length = 50;

    // Iterate through the alignments
    int num_alignments = 0;
    int num_sv_calls = 0;
    while (sam_read1(fp_in, bamHdr, bam1) >= 0) {

        // Get the tid and pos
        int32_t tid = bam1->core.tid;
        int32_t pos = bam1->core.pos;

        // Check if the data is sorted
        if (tid < last_tid || (tid == last_tid && pos < last_pos)) {
            sorted = false;

            // Stop checking if the data is sorted
            break;
        }
        last_tid = tid;
        last_pos = pos;

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
        }

        // Check if the insertion start and end positions have been set
        if (ins_start != -1 && ins_end != -1) {
            // Check if the insertion is longer than the minimum alignment length
            if (ins_end - ins_start > min_alignment_length) {

                // Add the insertion to the SV calls
                sv_calls[std::make_pair(chr, ins_start)] = std::make_pair(ins_end, 1);
                num_sv_calls++;
            }
        }

        // Check if the deletion start and end positions have been set
        if (del_start != -1 && del_end != -1) {
            // Check if the deletion is longer than the minimum alignment length
            if (del_end - del_start > min_alignment_length) {

                // Add the deletion to the SV calls
                sv_calls[std::make_pair(chr, del_start)] = std::make_pair(del_end, 2);
                num_sv_calls++;
            }
        }

        // Increment the number of alignments
        num_alignments++;

        // Print the progress as a wait cursor
        if (num_alignments % 100000 == 0) {
            std::cout << "\r" << num_alignments << " alignments processed";
            std::cout.flush();
        }
    }

    // Close the BAM file
    sam_close(fp_in);

    // Destroy the bam1_t
    bam_destroy1(bam1);

    // Destroy the bam_hdr_t
    bam_hdr_destroy(bamHdr);

    // Check if the data is sorted
    if (!sorted) {
        std::cerr << "ERROR: the data is not sorted" << std::endl;
        exit(1);
    }

    // Print the number of alignments processed
    std::cout << "\r" << num_alignments << " alignments processed" << std::endl;

    // Print the number of SV calls
    std::cout << num_sv_calls << " SV calls" << std::endl;

    // Return the SV calls
    return sv_calls;
}
