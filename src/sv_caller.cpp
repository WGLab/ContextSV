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


SVCaller::SVCaller(InputData& input_data)
{
    this->input_data = &input_data;
}

// Detect SVs and return SV type by start and end position
SVData SVCaller::run()
{
    // Get SV calls from split read alignments (primary and supplementary) and
    // directly from the CIGAR string
    SVData sv_calls = this->detectSVsFromSplitReads();
    //sv_calls.addSVCalls(split_read_calls);
    
    // Return the SV calls
    return sv_calls;
}

SVData SVCaller::detectSVsFromCIGAR(std::string chr, int32_t pos, uint32_t *cigar, int cigar_len)
{
    // Create a map of SV calls
    FASTAQuery ref_genome = this->input_data->getRefGenome();
    SVData sv_calls(ref_genome);

    // Track gaps between alignments for merging SVs
    int32_t merged_insertion_start = -1;
    int32_t last_insertion_end = -1;
    int32_t merged_deletion_start = -1;
    int32_t last_deletion_end = -1;

    // Store the read sequence
    std::string read_seq = "";

    // Get the end position
    int32_t end = pos;
    for (int i = 0; i < cigar_len; i++) {
        // Get the CIGAR operation
        int op = bam_cigar_op(cigar[i]);

        // Get the CIGAR operation length
        int op_len = bam_cigar_oplen(cigar[i]);

        // Process insertions
        if (op == BAM_CINS) {
            
            // Get the insertion start and end positions
            int32_t ins_start = end;
            int32_t ins_end = end + op_len;

            // Check if this is the first insertion
            if (merged_insertion_start == -1) {
                
                // Set the insertion start and end positions
                merged_insertion_start = ins_start;
                last_insertion_end = ins_end;
            
            // Check if the insertion is within the maximum distance from the
            // last insertion
            } else if (ins_start - last_insertion_end < this->max_indel_dist) {

                // Increment the last insertion's end position
                last_insertion_end = ins_end;

            } else {
                // Add the previous insertion to the SV calls
                //sv_calls.addSVCall(chr, merged_insertion_start, last_insertion_end, 0, "");


                // Print the SV type, position, and length if length is greater
                // than 20 kb
                if (last_insertion_end - merged_insertion_start > 20000) {
                    std::cout << "CIGAR INS\t" << merged_insertion_start << "\t" << last_insertion_end << "\t" << last_insertion_end - merged_insertion_start << std::endl;
                }
                //std::cout << "CIGAR INS\t" << merged_insertion_start << "\t" << last_insertion_end << "\t" << last_insertion_end - merged_insertion_start << std::endl;

                // Start a new insertion
                merged_insertion_start = ins_start;
                last_insertion_end = ins_end;

            }

        // Process deletions
        } else if (op == BAM_CDEL) {
                
            // Get the deletion start and end positions
            int32_t del_start = end;
            int32_t del_end = end + op_len;

            // Check if this is the first deletion
            if (merged_deletion_start == -1) {
                
                // Set the deletion start and end positions
                merged_deletion_start = del_start;
                last_deletion_end = del_end;
            
            // Check if the deletion is within the maximum distance from the
            // last deletion
            } else if (del_start - last_deletion_end < this->max_indel_dist) {

                // Increment the last deletion's end position
                last_deletion_end = del_end;

            } else {
                // Add the previous deletion to the SV calls
                //sv_calls.addSVCall(chr, merged_deletion_start, last_deletion_end, 1, "");
                
                // Print the SV type, position, and length if length is greater
                // than 20 kb
                if (last_deletion_end - merged_deletion_start > 20000) {
                    std::cout << "CIGAR DEL\t" << merged_deletion_start << "\t" << last_deletion_end << "\t" << last_deletion_end - merged_deletion_start << std::endl;
                }
                //std::cout << "CIGAR DEL\t" << merged_deletion_start << "\t" << last_deletion_end << "\t" << last_deletion_end - merged_deletion_start << std::endl;

                // Start a new deletion
                merged_deletion_start = del_start;
                last_deletion_end = del_end;

            }
        }

        // Increment the end position
        if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP || op == BAM_CEQUAL || op == BAM_CDIFF) {
            end += op_len;
        } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) {
            // Do nothing
        } else {
            std::cerr << "ERROR: Unknown CIGAR operation" << std::endl;
            exit(1);
        }
    }

    // Return the SV calls
    return sv_calls;
}

SVData SVCaller::detectSVsFromSplitReads()
{
    // Open the BAM file
    samFile *fp_in = sam_open(this->input_data->getBAMFilepath().c_str(), "r");
    if (fp_in == NULL) {
        std::cerr << "ERROR: failed to open " << this->input_data->getBAMFilepath() << std::endl;
        exit(1);
    }

    // Get the header
    bam_hdr_t *bamHdr = sam_hdr_read(fp_in);
    if (bamHdr == NULL) {
        std::cerr << "ERROR: failed to read header for " << this->input_data->getBAMFilepath() << std::endl;
        exit(1);
    }

    // Get the region
    std::string region = this->input_data->getRegion();

    // Get the index
    hts_idx_t *idx = sam_index_load(fp_in, this->input_data->getBAMFilepath().c_str());
    if (idx == NULL) {
        std::cerr << "ERROR: failed to load index for " << this->input_data->getBAMFilepath() << std::endl;
        exit(1);
    }

    // Set up the iterator
    hts_itr_t *itr = sam_itr_querys(idx, bamHdr, region.c_str());

    // Set up the read
    bam1_t *bam1 = bam_init1();

    // Create a map of primary and supplementary alignments by QNAME (query template name)
    int num_alignments = 0;
    int num_sv_calls = 0;
    FASTAQuery ref_genome = this->input_data->getRefGenome();
    SVData sv_calls(ref_genome);
    QueryMap primary_alignments;  // TODO: Add depth to primary alignments
    QueryMap supplementary_alignments;
    while (sam_itr_next(fp_in, itr, bam1) >= 0) {

        // Get the QNAME (query template name) for associating split reads
        std::string qname = bam_get_qname(bam1);

        // Skip if the QNAME is not equal to
        // A00739:176:H3NYTDSXY:3:1413:30273:31532
        // if (qname != "A00739:176:H3NYTDSXY:3:1413:30273:31532") {
        //     continue;
        // } else {
        //     std::cout << "Found Test QNAME" << std::endl;
        // }

        // Skip secondary and unmapped alignments
        if (bam1->core.flag & BAM_FSECONDARY || bam1->core.flag & BAM_FUNMAP) {
            // Do nothing

        // Skip alignments with low mapping quality
        } else if (bam1->core.qual < this->min_mapq) {
            // Do nothing
        } else {

            // Process primary alignments
            if (!(bam1->core.flag & BAM_FSUPPLEMENTARY)) {

                //std::cout << "Primary alignment" << std::endl;

                // Add the primary alignment to the map
                std::string chr = bamHdr->target_name[bam1->core.tid];
                int32_t start = bam1->core.pos;
                int32_t end = bam_endpos(bam1);
                int32_t depth = bam1->core.n_cigar;
                AlignmentData alignment(chr, start, end, depth);
                primary_alignments[qname].push_back(alignment);

                // Print the primary alignment position
                //std::cout << "Primary alignment" << std::endl;
                //std::cout << bamHdr->target_name[bam1->core.tid] << ":" << bam1->core.pos << "-" << bam_endpos(bam1) << std::endl;
            
                // Finally, call SVs directly from the CIGAR string
                //std::cout << "Calling SVs from CIGAR string..." << std::endl;                
                //SVData cigar_calls = this->detectSVsFromCIGAR(chr, bam1->core.pos, bam_get_cigar(bam1), bam1->core.n_cigar);
                //std::cout << "Complete." << std::endl;

                // Add the SV calls to the SV map
                //sv_calls.addSVCalls(cigar_calls);

            // Process supplementary alignments
            } else if (bam1->core.flag & BAM_FSUPPLEMENTARY) {

                std::cout << "Found supplementary alignment" << std::endl;

                // Add the supplementary alignment to the map
                std::string chr = bamHdr->target_name[bam1->core.tid];
                int32_t start = bam1->core.pos;
                int32_t end = bam_endpos(bam1);
                int32_t depth = bam1->core.n_cigar;
                AlignmentData alignment(chr, start, end, depth);

                // Add the supplementary alignment to the map
                supplementary_alignments[qname].push_back(alignment);

                // Print the supplementary alignment position
                //std::cout << "Supplementary alignment" << std::endl;
                //std::cout << bamHdr->target_name[bam1->core.tid] << ":" << bam1->core.pos << "-" << bam_endpos(bam1) << std::endl;
            }
        }

        // Increment the number of alignment records processed
        num_alignments++;
    }
    
    // Print the number of alignments processed
    std::cout << num_alignments << " alignments processed" << std::endl;

    // Loop through the map of primary alignments by QNAME and get the
    // supplementary alignments
    std::cout << "Running split read analysis..." << std::endl;
    for (const auto& entry : primary_alignments) {

        // Get the QNAME
        std::string qname = entry.first;

        // Get the first primary alignment
        //std::tuple<std::string, int, int> primary_alignment_tuple = entry.second[0];
        AlignmentData primary_alignment = entry.second[0];

        // Get the primary alignment chromosome
        std::string primary_chr = std::get<0>(primary_alignment);

        // Get the start and end positions of the primary alignment
        int32_t primary_start = std::get<1>(primary_alignment);
        int32_t primary_end = std::get<2>(primary_alignment);

        //std::cout << "Primary alignment at " << primary_chr << ":" << primary_start << "-" << primary_end << " with depth " << primary_depth << std::endl;

        // Identify potential duplications by looking for supplementary
        // alignments that overlap with the primary alignment


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
                std::cout << "Supplementary alignment on different chromosome" << std::endl;
                continue;
            }

            // Get the start and end positions of the supplementary alignment
            int32_t supp_start = std::get<1>(supp_alignment);
            int32_t supp_end = std::get<2>(supp_alignment);

            //std::cout << "Supplementary alignment at " << supp_chr << ":" << supp_start << "-" << supp_end << std::endl;

            // Determine the type of gap between alignments
            int32_t gap_start = 0;
            int32_t gap_end = 0;
            if (supp_start < primary_start && supp_end < primary_start) {
                // INDEL (Supplementary alignment is before primary
                // alignment with no overlap)
                // = Tandem duplication or translocation

                // // Print if within the region of interest 61135039-61911274
                // if (gap_start > 61000000 && gap_end < 62000000) {
                //     std::cout << "INDEL at " << supp_chr << ":" << gap_start << "-" << gap_end << std::endl;
                // } else {
                //     continue;
                // }

                // Get the gap start and end positions
                gap_start = supp_end;
                gap_end = primary_start;

                std::cout << "FWD GAP at " << supp_chr << ":" << gap_start << "-" << gap_end << std::endl;
                
            } else if (supp_start > primary_end && supp_end > primary_end) {
                // INDEL (Supplementary alignment is after primary alignment with no overlap)
                // = Deletion

                // // Print if within the region of interest 61135039-61911274
                // if (gap_start > 61000000 && gap_end < 62000000) {
                //     std::cout << "INDEL at " << supp_chr << ":" << gap_start << "-" << gap_end << std::endl;
                // } else {
                //     continue;
                // }

                // Get the gap start and end positions
                gap_start = primary_end;
                gap_end = supp_start;

                std::cout << "REV GAP at " << supp_chr << ":" << gap_start << "-" << gap_end << std::endl;

            } else {
                // OVERLAP (Supplementary alignment is within primary alignment or overlaps with primary alignment)
                // = Tandem duplication or translocation

                // // Print if within the region of interest 61135039-61911274
                // if (gap_start > 61000000 && gap_end < 62000000) {
                //     std::cout << "OVERLAP at " << supp_chr << ":" << gap_start << "-" << gap_end << std::endl;
                // } else {
                //     continue;
                // }

                // Get the gap start and end positions
                gap_start = std::min(primary_start, supp_start);
                gap_end = std::max(primary_end, supp_end);

                std::cout << "OVERLAP at " << supp_chr << ":" << gap_start << "-" << gap_end << std::endl;
            }

            // Check if the gap is larger than the minimum SV size
            if (gap_end - gap_start > this->min_sv_size) {

                // // Print if within the region of interest 61135039-61911274
                // if (gap_start > 61000000 && gap_end < 62000000) {
                //     std::cout << "INDEL at " << supp_chr << ":" << gap_start << "-" << gap_end << std::endl;
                // } else {
                //     continue;
                // }

                // Add the gap to the SV calls
                // Set the SV type to -1 for now, will be labeled later using CNV calls
                sv_calls.addSVCall(supp_chr, gap_start, gap_end, -1, ".");
                
                //std::cout << "Added SV call at " << supp_chr << ":" << gap_start << "-" << gap_end << std::endl;
                std::cout << "SV size = " << gap_end - gap_start << std::endl;

                // Increment the number of SV calls
                num_sv_calls++;
            }
        }
    }

    // Destroy objects and close files
    bam_destroy1(bam1);  // Destroy the read 
    hts_itr_destroy(itr);  // Destroy the iterator
    sam_close(fp_in);  // Close the BAM file
    bam_hdr_destroy(bamHdr);  // Destroy the header
    hts_idx_destroy(idx);  // Destroy the index

    // Print the number of SV calls
    std::cout << num_sv_calls << " SV calls" << std::endl;

    // Return the SV calls
    return sv_calls;
}
