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
    // SVMap sv_calls = this->detectSVsFromCIGAR();

    // Get SV calls from split read alignments (primary and supplementary)
    SVMap sv_calls = this->detectSVsFromSplitReads();

    // // Merge SV calls
    // SVMap merged_sv_calls = this->mergeSVs(sv_calls);

    // // Get the SV type for each SV call
    // SVMap sv_calls_with_type = this->predictSVType(merged_sv_calls,
    // cnv_calls);
    
    // Return the SV calls
    return sv_calls;
}

SVMap SVCaller::detectSVsFromCIGAR()
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
        
        // Skip alignments that are not primary or supplementary
        if (!(bam1->core.flag & BAM_FSECONDARY) && !(bam1->core.flag & BAM_FSUPPLEMENTARY)) {
            continue;
        }

        // Skip alignments with mapping quality less than min_mapq
        if (bam1->core.qual < min_mapq) {
            continue;
        }

        // Look for small insertions and deletions in the CIGAR string for
        // primary and supplementary alignments
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

                        // Print the deletion start and end positions and size
                        // if within the region chr6:15482095-15729539
                        // Tolerate 1 kb difference in start and end positions
                        int32_t region_start = 15482095;
                        int32_t region_end = 15729539;
                        if (chr == "chr6" && merged_del_start >= region_start - 0 && merged_del_end <= region_end + 0) {
                            std::cout << "Deletion at " << chr << ":" << merged_del_start << "-" << merged_del_end << std::endl;
                            std::cout << "Deletion size = " << merged_del_end - merged_del_start << std::endl;
                        }
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

SVMap SVCaller::detectSVsFromSplitReads()
{
    // Open the BAM file
    samFile *fp_in = sam_open(this->common.getBAMFilepath().c_str(), "r");
    if (fp_in == NULL) {
        std::cerr << "ERROR: failed to open " << this->common.getBAMFilepath() << std::endl;
        exit(1);
    }

    // Get the header
    bam_hdr_t *bamHdr = sam_hdr_read(fp_in);
    if (bamHdr == NULL) {
        std::cerr << "ERROR: failed to read header for " << this->common.getBAMFilepath() << std::endl;
        exit(1);
    }

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

    // Iterate through the alignments in the region and get the gaps between
    // split reads
    
    SVMap sv_calls;
    int num_alignments = 0;
    int num_sv_calls = 0;

    // Create a map of primary and supplementary alignments by QNAME and value
    // of chr:start-end
    // Note: Multiple primary alignments can have the same QNAME when alignment
    // quality is similar ex. in repetitive regions
    QueryMap primary_alignments;
    QueryMap supplementary_alignments;
    while (sam_itr_next(fp_in, itr, bam1) >= 0) {

        // Get the QNAME (query template name) for associating split reads
        std::string qname = bam_get_qname(bam1);

        // Skip if the QNAME is not equal to
        // A00739:176:H3NYTDSXY:3:1413:30273:31532
        if (qname != "A00739:176:H3NYTDSXY:3:1413:30273:31532") {
            continue;
        } else {
            std::cout << "Found Test QNAME" << std::endl;
        }

        // Check that the alignment is primary
        if (!(bam1->core.flag & BAM_FSECONDARY) && !(bam1->core.flag & BAM_FSUPPLEMENTARY)) {
            // Add the primary alignment to the map
            std::string chr = bamHdr->target_name[bam1->core.tid];
            int32_t start = bam1->core.pos;
            int32_t end = bam_endpos(bam1);
            
            // Add the primary alignment to the map
            AlignmentLocation location(chr, start, end);

            // Add the primary alignment to the map
            primary_alignments[qname].push_back(location);

            // Add the primary alignment to the map
            //query_map.addAlignment(chr, start, end, 0);

            // Print the primary alignment position
            std::cout << "Primary alignment" << std::endl;
            std::cout << bamHdr->target_name[bam1->core.tid] << ":" << bam1->core.pos << "-" << bam_endpos(bam1) << std::endl;

            // // Print the primary alignment position if equal to
            // // chr6:15,482,028-15,482,097
            // if (bamHdr->target_name[bam1->core.tid] == "chr6" && bam1->core.pos == 15482028) {
            //     std::cout << "Primary alignment" << std::endl;
            //     std::cout << bamHdr->target_name[bam1->core.tid] << ":" << bam1->core.pos << "-" << bam_endpos(bam1) << std::endl;
            
            //     // Print the QNAME
            //     std::cout << qname << std::endl;
            // }

            // // Get the supplementary alignment positions
            // std::vector<bam1_t*> supp_alignments = supp_reads[qname];
            // for (const auto& supp_alignment : supp_alignments) {
            //     std::cout << "Supplementary alignment" << std::endl;
            //     std::cout << bamHdr->target_name[supp_alignment->core.tid] << ":" << supp_alignment->core.pos << "-" << bam_endpos(supp_alignment) << std::endl;
            // }

        }

        // Check that the alignment is supplementary
        else if (bam1->core.flag & BAM_FSUPPLEMENTARY) {
            // Add the supplementary alignment to the map
            std::string chr = bamHdr->target_name[bam1->core.tid];
            int32_t start = bam1->core.pos;
            int32_t end = bam_endpos(bam1);
            //alignments.addAlignment(chr, start, end, 1);

            // Add the supplementary alignment to the map
            AlignmentLocation location(chr, start, end);

            // Add the supplementary alignment to the map
            supplementary_alignments[qname].push_back(location);

            // Print the supplementary alignment position
            std::cout << "Supplementary alignment" << std::endl;
            std::cout << bamHdr->target_name[bam1->core.tid] << ":" << bam1->core.pos << "-" << bam_endpos(bam1) << std::endl;

            // // Print the supplementary alignment position if equal to
            // // chr6:15,729,539-15,729,622
            // if (bamHdr->target_name[bam1->core.tid] == "chr6" && bam1->core.pos == 15729539) {
            //     std::cout << "Supplementary alignment" << std::endl;
            //     std::cout << bamHdr->target_name[bam1->core.tid] << ":" << bam1->core.pos << "-" << bam_endpos(bam1) << std::endl;

            //     // Print the QNAME
            //     std::cout << qname << std::endl;
            // }
        }

        // Increment the number of alignment records processed
        num_alignments++;

    }
    
    // Print the number of alignments processed
    std::cout << num_alignments << " alignment records processed" << std::endl;

    // Iterate through the primary alignments and get the gaps between split
    // reads
    std::cout << "Calling SVs from split reads..." << std::endl;

    // Loop through the map of primary alignments by QNAME
    //std::map<std::string, std::vector<std::tuple<std::string, int, int>>> primary_reads = alignments.getPrimaryAlignments();
    for (const auto& entry : primary_alignments) {

        // Get the QNAME
        std::string qname = entry.first;

        // Get the first primary alignment
        //std::tuple<std::string, int, int> primary_alignment_tuple = entry.second[0];
        AlignmentLocation primary_alignment = entry.second[0];

        // Get the primary alignment chromosome
        std::string primary_chr = std::get<0>(primary_alignment);

        // Get the start and end positions of the primary alignment
        int32_t primary_start = std::get<1>(primary_alignment);
        int32_t primary_end = std::get<2>(primary_alignment);

        std::cout << "Primary alignment at " << primary_chr << ":" << primary_start << "-" << primary_end << std::endl;

        // Get the supplementary alignments
        //std::vector<std::tuple<std::string, int, int>> supp_alignments = alignments.getSupplementaryAlignments(qname);
        AlignmentVector supp_alignments = supplementary_alignments[qname];

        // Get the gaps between split reads
        // Loop through the supplementary alignments
        for (const auto& supp_alignment : supp_alignments) {

            // Get the supplementary alignment chromosome
            std::string supp_chr = std::get<0>(supp_alignment);

            // Skip supplementary alignments that are on a different chromosome
            // for now (TODO: Use for translocations)
            if (primary_chr != supp_chr) {
                continue;
            }

            // Get the start and end positions of the supplementary alignment
            int32_t supp_start = std::get<1>(supp_alignment);
            int32_t supp_end = std::get<2>(supp_alignment);

            std::cout << "Supplementary alignment at " << supp_chr << ":" << supp_start << "-" << supp_end << std::endl;

            // Determine if it is an insertion or deletion (Insertion = overlap
            // between primary and supplementary alignments, Deletion = gap
            // between primary and supplementary alignments)
            int gap_type = 0;  // 0 = insertion, 1 = deletion
            int32_t gap_start = 0;
            int32_t gap_end = 0;
            if (supp_start < primary_start && supp_end < primary_start) {
                // Insertion (Supplementary alignment is before primary
                // alignment with no overlap)
                gap_type = 0;

                // Get the gap start and end positions
                gap_start = supp_end;
                gap_end = primary_start;

                std::cout << "Insertion at " << supp_chr << ":" << gap_start << "-" << gap_end << std::endl;
                
            } else if (supp_start > primary_end && supp_end > primary_end) {
                // Insertion (Supplementary alignment is after primary alignment with no overlap)
                gap_type = 0;

                // Get the gap start and end positions
                gap_start = primary_end;
                gap_end = supp_start;

                std::cout << "Insertion at " << supp_chr << ":" << gap_start << "-" << gap_end << std::endl;

            } else {
                // Deletion (Supplementary alignment is within primary alignment or overlaps with primary alignment)
                gap_type = 1;

                // Get the gap start and end positions
                gap_start = std::min(primary_start, supp_start);
                gap_end = std::max(primary_end, supp_end);

                std::cout << "Deletion at " << supp_chr << ":" << gap_start << "-" << gap_end << std::endl;
            }

            // Check if the gap is larger than the minimum SV size
            if (gap_end - gap_start > this->min_sv_size) {

                // Add the gap to the SV calls
                sv_calls.addSVCall(supp_chr, gap_start, gap_end, gap_type);

                //std::cout << "Added SV call at " << supp_chr << ":" << gap_start << "-" << gap_end << std::endl;
                std::cout << "SV size = " << gap_end - gap_start << std::endl;

                // Increment the number of SV calls
                num_sv_calls++;
            }
        }

        // // Destroy the primary alignment
        // bam_destroy1(primary_alignment);

        // // Destroy the supplementary alignments
        // for (const auto& supp_alignment : supp_alignments) {
        //     bam_destroy1(supp_alignment);
        // }
    }

    // Destroy the alignment
    bam_destroy1(bam1);

    // Destroy the iterator
    hts_itr_destroy(itr);

    // Close the BAM file
    sam_close(fp_in);

    // Destroy the header
    bam_hdr_destroy(bamHdr);

    // Destroy the index
    hts_idx_destroy(idx);

    // Print the number of SV calls
    std::cout << num_sv_calls << " SV calls" << std::endl;

    // Return the SV calls
    return sv_calls;
}
