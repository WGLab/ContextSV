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
/// @endcond

SVCaller::SVCaller(InputData& input_data)
{
    this->input_data = &input_data;
}

// Main function for SV detection
SVData SVCaller::run()
{
    // Get SV calls from split read alignments (primary and supplementary) and
    // directly from the CIGAR string
    SVData sv_calls = this->detectSVsFromSplitReads();
    
    // Return the SV calls
    return sv_calls;
}

// Detect SVs from the CIGAR string of a read alignment.
void SVCaller::detectSVsFromCIGAR(SVData& sv_calls, std::string chr, int32_t pos, uint32_t *cigar, int cigar_len, bool debug_mode)
{
    // Create a map of SV calls
    FASTAQuery ref_genome = this->input_data->getRefGenome();

    // Store the read sequence
    std::string read_seq = "";

    // Print the start position
    if (debug_mode) {
        std::cout << "CIGAR start position: " << chr << ":" << pos << std::endl;
    }

    // Loop through the CIGAR string and detect insertions and deletions in
    // reference coordinates
        //pos += 1;  // Update the position +=1 for 1-based coordinates
    // For BAM files, the reference starts at 0. For SAM, 1. POS is the leftmost position of where the alignment maps to the reference:
    // https://genome.sph.umich.edu/wiki/SAM
    for (int i = 0; i < cigar_len; i++) {

        // Get the CIGAR operation
        int op = bam_cigar_op(cigar[i]);

        // Get the CIGAR operation length
        int op_len = bam_cigar_oplen(cigar[i]);
        
        // Check if the CIGAR operation is an insertion
        if (op == BAM_CINS) {

            // Add the SV if greater than the minimum SV size
            if (op_len >= this->min_sv_size) {
                // Add the insertion to the SV calls
                //sv_calls.add(chr, pos, pos + op_len, 0, ".", "CIGARINS");
                sv_calls.add(chr, pos, pos + op_len, 3, ".", "CIGARINS");

                if (debug_mode) {
                    std::cout << "Found CIGAR insertion, LEN=" << op_len << " at " << chr << ":" << pos << "-" << pos+op_len << std::endl;
                }
            }

        // Check if the CIGAR operation is a deletion
        } else if (op == BAM_CDEL) {

            // Add the SV if greater than the minimum SV size
            if (op_len >= this->min_sv_size) {
                
                // Add the deletion to the SV calls
                sv_calls.add(chr, pos, pos + op_len, 0, ".", "CIGARDEL");

                // // Print if 21:37857255-37857617
                // if (chr == "21" && pos == 37857255 && pos + op_len == 37857617) {
                //     std::cout << "[CIGAR FIND] Added 21:37857255-37857617 type: " << 0 << std::endl;
                // }

                // // Print if 21:42911938-42912387
                // if (chr == "21" && pos == 42911938 && pos + op_len == 42912387) {
                //     std::cout << "[CIGAR FIND] Added 21:42911938-42912387 type: " << 0 << std::endl;
                // }

                if (debug_mode) {
                    std::cout << "Found CIGAR deletion, LEN=" << op_len << " at " << chr << ":" << pos << "-" << pos+op_len << std::endl;
                }
            }
        }

        // Update the reference coordinate based on the CIGAR operation (M, D,
        // N, =, X)
        // https://samtools.github.io/hts-specs/SAMv1.pdf
        if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP || op == BAM_CEQUAL || op == BAM_CDIFF) {
            pos += op_len;
        } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP || op == BAM_CPAD) {
            // Do nothing
        } else {
            std::cerr << "ERROR: Unknown CIGAR operation " << op << std::endl;
            exit(1);
        }
    }
}

// Detect SVs from split read alignments (primary and supplementary) and
// directly from the CIGAR string
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
    //int num_sv_calls = 0;
    FASTAQuery ref_genome = this->input_data->getRefGenome();
    SVData sv_calls(ref_genome);
    QueryMap primary_alignments;  // TODO: Add depth to primary alignments
    QueryMap supplementary_alignments;
    std::cout << "Reading alignments..." << std::endl;
    while (sam_itr_next(fp_in, itr, bam1) >= 0) {

        // Get the QNAME (query template name) for associating split reads
        std::string qname = bam_get_qname(bam1);

        // If the read name is equal to 87c5f638-7716-491f-bee2-d04110f6e743
        // then print CIGAR debug information
        bool debug_cigar = false;
        //if (qname == "87c5f638-7716-491f-bee2-d04110f6e743") {
        // if (qname == "5bd073fe-e4b0-4abc-9535-229f13249578") {
        //    debug_cigar = true;
        //    std::cout << "Found Test QNAME " << qname << std::endl;
        //} else {
        //    // Skip for debugging purposes
        //    continue;
        //}

        //Skip if the QNAME is not equal to
        //A00739:176:H3NYTDSXY:3:1413:30273:31532
        // if (qname == "08033753-0fa2-4321-a300-939df26c5af1" || qname == "9746d47c-75ee-40e1-ba39-37ff5dc6e028") {
        //     std::cout << "Found Test QNAME" << std::endl;
        // } else {
        //     //std::cout << "Skipping QNAME " << qname << std::endl;
        //     continue;
        // }

        // Skip secondary and unmapped alignments
        if (bam1->core.flag & BAM_FSECONDARY || bam1->core.flag & BAM_FUNMAP) {
            // Do nothing

        // Skip alignments with low mapping quality
        } else if (bam1->core.qual < this->min_mapq) {
            // Do nothing
            // std::cout << "Skipping alignment with low mapping quality" << std::endl;
            // std::cout << bam1->core.qual << " < " << this->min_mapq << std::endl;
        } else {

            // Process primary alignments
            if (!(bam1->core.flag & BAM_FSUPPLEMENTARY)) {

                // Add the primary alignment to the map
                std::string chr = bamHdr->target_name[bam1->core.tid];
                int64_t start = bam1->core.pos;
                int64_t end = bam_endpos(bam1);
                int32_t cigar_len = bam1->core.n_cigar;
                int depth = 0;  // Placeholder for now
                AlignmentData alignment(chr, start, end, depth);
                primary_alignments[qname].push_back(alignment);

                // Print the entire CIGAR string
                if (debug_cigar) {
                    std::cout << "FULL CIGAR: " << std::endl;
                    for (uint32_t i = 0; i < bam1->core.n_cigar; i++) {
                        std::cout << bam1->core.n_cigar << bam_cigar_opchr(bam1->core.n_cigar) << bam_cigar_oplen(bam1->core.n_cigar) << " ";
                    }
                    std::cout << std::endl;
                }
            
                // Finally, call SVs directly from the CIGAR string
                //int prev_sv_count = sv_calls.size();
                if (!this->input_data->getDisableCIGAR()) {

                    //std::cout << "Calling SVs from CIGAR string" << std::endl;
                    this->detectSVsFromCIGAR(sv_calls, chr, start, bam_get_cigar(bam1), cigar_len, debug_cigar);
                }
                //this->detectSVsFromCIGAR(sv_calls, chr, bam1->core.pos, bam_get_cigar(bam1), bam1->core.n_cigar);
                //int curr_sv_count = sv_calls.size();
                // if (curr_sv_count > prev_sv_count) {
                //     std::cout << "Found " << curr_sv_count - prev_sv_count << "CIGAR SVs" << std::endl;
                // }

            // Process supplementary alignments
            } else if (bam1->core.flag & BAM_FSUPPLEMENTARY) {

                //std::cout << "Found supplementary alignment" << std::endl;

                // Add the supplementary alignment to the map
                std::string chr = bamHdr->target_name[bam1->core.tid];
                int32_t start = bam1->core.pos;
                int32_t end = bam_endpos(bam1);
                int32_t depth = bam1->core.n_cigar;
                AlignmentData alignment(chr, start, end, depth);

                // Add the supplementary alignment to the map
                supplementary_alignments[qname].push_back(alignment);
            }
        }

        // Increment the number of alignment records processed
        num_alignments++;

        // Print the number of alignments processed every 10 thousand
        if (num_alignments % 10000 == 0) {
            std::cout << num_alignments << " alignments processed" << std::endl;
        }
    }
    
    // Print the number of alignments processed
    std::cout << num_alignments << " alignments processed" << std::endl;

    // Loop through the map of primary alignments by QNAME and find gaps and
    // overlaps from supplementary alignments
    std::cout << "Running split read analysis..." << std::endl;
    for (const auto& entry : primary_alignments) {

        // Get the QNAME
        std::string qname = entry.first;

        // Get the first primary alignment
        AlignmentData primary_alignment = entry.second[0];

        // Get the primary alignment chromosome
        std::string primary_chr = std::get<0>(primary_alignment);

        // Get the start and end positions of the primary alignment
        int32_t primary_start = std::get<1>(primary_alignment);
        int32_t primary_end = std::get<2>(primary_alignment);

        // Loop through the supplementary alignments and find gaps and overlaps
        AlignmentVector supp_alignments = supplementary_alignments[qname];
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


            // Determine the type of gap between alignments
            if (supp_start < primary_start && supp_end < primary_start) {
                // Gap with supplementary before primary:
                // [supp_start] [supp_end] -- [primary_start] [primary_end]

                // Use the gap ends as the SV endpoints (Deletion)
                if (primary_start - supp_end >= this->min_sv_size) {
                    sv_calls.add(supp_chr, supp_end, primary_start, -1, ".", "GAPINNER_1");
                }

                // Use the alignment ends as the SV endpoints (Insertion)
                if (primary_end - supp_start >= this->min_sv_size) {
                    sv_calls.add(supp_chr, supp_start, primary_end, -1, ".", "GAPOUTER_1");
                }

                //std::cout << "FWD GAP at " << supp_chr << ":" << gap_start << "-" << gap_end << std::endl;
                
            } else if (supp_start > primary_end && supp_end > primary_end) {
                // Gap with supplementary after primary:
                // [primary_start] [primary_end] -- [supp_start] [supp_end]

                // Use the gap ends as the SV endpoints (Deletion)
                if (supp_start - primary_end >= this->min_sv_size) {
                    sv_calls.add(supp_chr, primary_end, supp_start, -1, ".", "GAPINNER_2");
                }

                // Use the alignment ends as the SV endpoints (Insertion)
                if (supp_end - primary_start >= this->min_sv_size) {
                    sv_calls.add(supp_chr, primary_start, supp_end, -1, ".", "GAPOUTER_2");
                }

                //std::cout << "REV GAP at " << supp_chr << ":" << gap_start << "-" << gap_end << std::endl;

            } else {
                // Overlap between alignments:
                // [supp_start] [primary_start] [supp_end] [primary_end]
                // [primary_start] [supp_start] [primary_end] [supp_end]

                // Use the overlap ends as the SV endpoints (Deletion)
                if (supp_end - primary_start >= this->min_sv_size) {
                    sv_calls.add(supp_chr, primary_start, supp_end, -1, ".", "OVERLAP_1");
                } else if (primary_end - supp_start > this->min_sv_size) {
                    sv_calls.add(supp_chr, supp_start, primary_end, -1, ".", "OVERLAP_2");
                }

                // Use the alignment ends as the SV endpoints (Insertion)
                if (primary_end - supp_start >= this->min_sv_size) {
                    sv_calls.add(supp_chr, supp_start, primary_end, -1, ".", "OVERLAP_3");
                } else if (supp_end - primary_start > this->min_sv_size) {
                    sv_calls.add(supp_chr, primary_start, supp_end, -1, ".", "OVERLAP_4");
                }

                //std::cout << "OVERLAP at " << supp_chr << ":" << gap_start << "-" << gap_end << std::endl;
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
    std::cout << sv_calls.size() << " SV calls" << std::endl;

    // Return the SV calls
    return sv_calls;
}
