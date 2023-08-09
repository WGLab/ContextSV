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

    // Get SV calls from long read alignments
    SVMap sv_calls = this->detectSVs();

    // // Merge SV calls
    // SVMap merged_sv_calls = this->mergeSVs(sv_calls);

    // // Get the SV type for each SV call
    // SVMap sv_calls_with_type = this->predictSVType(merged_sv_calls,
    // cnv_calls);
    
    // Return the SV calls
    return sv_calls;
}

SVMap SVCaller::detectSVs()
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

    // Iterate through the alignments in the region and get the merged SV calls
    SVMap sv_calls;
    int num_alignments = 0;
    int num_sv_calls = 0;

    // Track the distance between consecutive insertions or deletions
    // This is used to merge consecutive insertions or deletions into a single
    // SV by setting the maximum distance between consecutive insertions or
    // deletions.
    int32_t merged_del_start = -1;  // Start position of the current deletion to be merged
    int32_t merged_del_end = -1;  // End position of the current deletion to be merged

    // TODO: Insertion distance

    while (sam_itr_next(fp_in, itr, bam1) >= 0) {
        
        // Skip secondary alignments
        if (bam1->core.flag & BAM_FSECONDARY) {            
            continue;
        }

        // Skip alignments with mapping quality less than min_mapq
        if (bam1->core.qual < min_mapq) {
            continue;
        }

        // Look for insertions and deletions in the CIGAR string
        int32_t ref_pos = bam1->core.pos;  // Reference coordinate
        uint32_t *cigar = bam_get_cigar(bam1);
        std::string chr = bamHdr->target_name[bam1->core.tid];
        for (uint32_t i = 0; i < bam1->core.n_cigar; i++) {

            // Get the CIGAR operation
            int cigar_op = bam_cigar_op(cigar[i]);

            // Get the CIGAR operation length
            int cigar_op_len = bam_cigar_oplen(cigar[i]);

            // Check if the CIGAR operation is an insertion
            if (cigar_op == BAM_CINS) {
                // TODO
            }

            // Check if the CIGAR operation is a deletion
            else if (cigar_op == BAM_CDEL) {

                // Get the deletion start and end positions
                int32_t del_start = ref_pos;
                int32_t del_end = ref_pos + cigar_op_len;

                // Check if this is the first deletion
                if (merged_del_start == -1) {

                    // Set the merged deletion start and end positions
                    merged_del_start = del_start;
                    merged_del_end = del_end;

                // Check if the distance between the current deletion and the
                // previous merged deletion is less than the maximum distance
                // between consecutive deletions. If so, merge the current
                // deletion with the previous merged deletion.
                // Note: this assumes that the CIGAR string is sorted by
                // reference position.
                // Also note: If there is overlap between the current deletion
                // and the previous merged deletion, the distance between the
                // two deletions will also be less than the maximum distance
                // between consecutive deletions.
                } else if (del_start - merged_del_end <= this->max_indel_dist) {

                    // Update the merged deletion end position
                    merged_del_end = del_end;
                
                // If the deletion is not within the maximum distance between
                // consecutive deletions, save the merged deletion and start a
                // new merged deletion
                } else {

                    // Check if the merged deletion is longer than the minimum SV size
                    if (merged_del_end - merged_del_start > this->min_sv_size) {

                        // Add the merged deletion to the SV calls
                        sv_calls.addSVCall(chr, merged_del_start, merged_del_end, 2);
                        num_sv_calls++;

                        // std::cout << "Added deletion at " << chr << ":" << merged_del_start << "-" << merged_del_end << std::endl;
                        // std::cout << "Deletion size = " << merged_del_end - merged_del_start << std::endl;

                        // // Print deletion size if larger than 100 kb
                        // if (merged_del_end - merged_del_start > 1000) {
                        //     std::cout << "Large deletion size = " << merged_del_end - merged_del_start << std::endl;
                        //     std::cout << "Deletion at " << chr << ":" << merged_del_start << "-" << merged_del_end << std::endl;
                        // }
                    }

                    // Reset the merged deletion start and end positions
                    merged_del_start = del_start;
                    merged_del_end = del_end;

                }
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
        }

        // Increment the number of alignment records processed
        num_alignments++;

        // // Print the progress as a wait cursor
        // if (num_alignments % 100000 == 0) {
        //     std::cout << "\r" << num_alignments << " alignment records processed";
        //     std::cout.flush();
        // }
    }

    // Print the number of alignments processed
    std::cout << num_alignments << " alignment records processed" << std::endl;

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

SVMap SVCaller::mergeSVs(SVMap sv_calls)
{
    // TODO
    return sv_calls;
}

SVMap SVCaller::predictSVType(SVMap sv_calls, CNVMap cnv_calls)
{
    return SVMap();
}
