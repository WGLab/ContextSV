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
/// @endcond

int SVCaller::readNextAlignment(samFile *fp_in, hts_itr_t *itr, bam1_t *bam1, std::mutex &mtx_bam)
{
    // Lock the mutex while reading the next alignment
    std::lock_guard<std::mutex> lock(mtx_bam);

    // Read the next alignment
    int ret = sam_itr_next(fp_in, itr, bam1);

    // Return the result of reading the next alignment
    return ret;
}

void SVCaller::detectSVsFromRegion(std::string region, SVData &sv_calls, samFile *fp_in, bam_hdr_t *bamHdr, hts_idx_t *idx, std::mutex &mtx_sv_calls, std::mutex &mtx_bam)
{
    // Create a read and iterator for the region
    bam1_t *bam1 = bam_init1();
    hts_itr_t *itr = sam_itr_querys(idx, bamHdr, region.c_str());

    // Loop through the alignments
    // Create a map of primary and supplementary alignments by QNAME (query template name)
    int num_alignments = 0;
    QueryMap primary_alignments;  // TODO: Add depth to primary alignments
    QueryMap supplementary_alignments;
    std::cout << "Reading alignments..." << std::endl;
    while (readNextAlignment(fp_in, itr, bam1, mtx_bam) >= 0) {

        // Skip secondary and unmapped alignments
        if (bam1->core.flag & BAM_FSECONDARY || bam1->core.flag & BAM_FUNMAP) {
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

                AlignmentData alignment(chr, start, end, ".");
                primary_alignments[qname].push_back(alignment);
            
                // Call SVs directly from the CIGAR string
                if (!this->input_data->getDisableCIGAR()) {

                    //std::cout << "Calling SVs from CIGAR string" << std::endl;
                    this->detectSVsFromCIGAR(bamHdr, bam1, sv_calls, mtx_sv_calls);
                }

            // Process supplementary alignments
            } else if (bam1->core.flag & BAM_FSUPPLEMENTARY) {

                // Add the supplementary alignment to the map
                std::string chr = bamHdr->target_name[bam1->core.tid];
                int32_t start = bam1->core.pos;
                int32_t end = bam_endpos(bam1);
                AlignmentData alignment(chr, start, end, ".");

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
            // for now (TODO: Use for identifying trans-chromosomal SVs such as
            // translocations)
            if (primary_chr != supp_chr) {
                std::cout << "Supplementary alignment on different chromosome" << std::endl;
                continue;
            }

            // Get the start and end positions of the supplementary alignment
            int32_t supp_start = std::get<1>(supp_alignment);
            int32_t supp_end = std::get<2>(supp_alignment);

            // Gap analysis (deletion or duplication)
            if (supp_start < primary_start && supp_end < primary_start) {
                // Gap with supplementary before primary:
                // [supp_start] [supp_end] -- [primary_start] [primary_end]

                // Use the gap ends as the SV endpoints (Possible insertion
                // (dup) or deletion)
                if (primary_start - supp_end >= this->min_sv_size) {

                    // Likely deletion
                    // SNP CNV predictions will be used to classify the SV type.
                    sv_calls.add(supp_chr, supp_end+1, primary_start+1, SVData::UNKNOWN, ".", "GAPINNER_A", mtx_sv_calls);
                }

                // [TEST] Use the alignment ends as the SV endpoints
                if (primary_end - supp_start >= this->min_sv_size) {
                    sv_calls.add(supp_chr, supp_start+1, primary_end+1, SVData::UNKNOWN, ".", "GAPOUTER_A", mtx_sv_calls);
                }

                
            } else if (supp_start > primary_end && supp_end > primary_end) {
                // Gap with supplementary after primary:
                // [primary_start] [primary_end] -- [supp_start] [supp_end]

                // Use the gap ends as the SV endpoints (Possible deletion)
                if (supp_start - primary_end >= this->min_sv_size) {
                    // // Add duplication
                    // sv_calls.add(supp_chr, primary_end+1, supp_start+1, SV_TYPE_DUP, ".", "GAPINNER_2", mtx_sv_calls);

                    // Likely deletion
                    // SNP CNV predictions will be used to classify the SV type.
                    sv_calls.add(supp_chr, primary_end+1, supp_start+1, SVData::UNKNOWN, ".", "GAPINNER_B", mtx_sv_calls);
                }

                // [TEST] Use the alignment ends as the SV endpoints
                if (supp_end - primary_start >= this->min_sv_size) {
                    sv_calls.add(supp_chr, primary_start+1, supp_end+1, SVData::UNKNOWN, ".", "GAPOUTER_B", mtx_sv_calls);
                }
            }
        }
    }

    // Destroy the iterator
    hts_itr_destroy(itr);

    // Destroy the read
    bam_destroy1(bam1);

    // Print the number of SV calls
    std::cout << sv_calls.size() << " SV calls from " << region << std::endl;
}

SVCaller::SVCaller(InputData &input_data)
{
    this->input_data = &input_data;
}

// Main function for SV detection
void SVCaller::run(SVData& sv_calls)
{
    // Get SV calls from split read alignments (primary and supplementary) and
    // directly from the CIGAR string
    this->detectSVsFromSplitReads(sv_calls);
}

// Detect SVs from the CIGAR string of a read alignment.
void SVCaller::detectSVsFromCIGAR(bam_hdr_t* header, bam1_t* alignment, SVData& sv_calls, std::mutex& mtx)
{
    // Get the chromosome
    std::string chr = header->target_name[alignment->core.tid];

    // Get the position
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
    int32_t ref_pos;
    int32_t ref_end;
    for (int i = 0; i < cigar_len; i++) {

        // Get the CIGAR operation
        int op = bam_cigar_op(cigar[i]);

        // Get the CIGAR operation length
        int op_len = bam_cigar_oplen(cigar[i]);
        
        // Check if the CIGAR operation is an insertion
        if (op == BAM_CINS) {

            // Add the SV if greater than the minimum SV size
            if (op_len >= this->min_sv_size) {

                // Get the sequence of the insertion from the query
                std::string ins_seq_str = "";
                uint8_t* seq_ptr = bam_get_seq(alignment);
                for (int j = 0; j < op_len; j++) {
                    ins_seq_str += seq_nt16_str[bam_seqi(seq_ptr, query_pos + j)];
                }

                // Add the insertion to the SV calls (1-based)
                //sv_calls.add(chr, pos, pos + op_len, 3, ins_seq_str,
                //"CIGARINS");
                ref_pos = pos+1;
                ref_end = ref_pos + op_len -1;
                sv_calls.add(chr, ref_pos, ref_end, SVData::INS, ins_seq_str, "CIGARINS", mtx);
            }

        // Check if the CIGAR operation is a deletion
        } else if (op == BAM_CDEL) {

            // Add the SV if greater than the minimum SV size
            if (op_len >= this->min_sv_size) {
                
                // Add the deletion to the SV calls (1-based)
                ref_pos = pos+1;
                ref_end = ref_pos + op_len -1;
                sv_calls.add(chr, ref_pos, ref_end, SVData::DEL, ".", "CIGARDEL", mtx);
            }

        // Check if the CIGAR operation is a soft clip
        } else if (op == BAM_CSOFT_CLIP) {

            // Update the clipped base support
            sv_calls.updateClippedBaseSupport(chr, pos);
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
}

// Detect SVs from split read alignments (primary and supplementary) and
// directly from the CIGAR string
SVData SVCaller::detectSVsFromSplitReads(SVData& sv_calls)
{
    // Open the BAM file
    std::string bam_filepath = this->input_data->getLongReadBam();
    samFile *fp_in = sam_open(bam_filepath.c_str(), "r");
    if (fp_in == NULL) {
        std::cerr << "ERROR: failed to open " << bam_filepath << std::endl;
        exit(1);
    }

    // Get the header
    bam_hdr_t *bamHdr = sam_hdr_read(fp_in);
    if (bamHdr == NULL) {
        std::cerr << "ERROR: failed to read header for " << bam_filepath << std::endl;
        exit(1);
    }

    // Get the region data
    bool whole_genome = this->input_data->getWholeGenome();
    std::vector<std::string> regions;
    if (whole_genome) {
        regions = this->input_data->getRefGenome().getChromosomes();
    } else {
        regions.push_back(this->input_data->getRegion());
    }

    // Get the index
    hts_idx_t *idx = sam_index_load(fp_in, bam_filepath.c_str());
    if (idx == NULL) {
        std::cerr << "ERROR: failed to load index for " << bam_filepath << std::endl;
        exit(1);
    }

    // Loop through the regions, creating a thread for each region
    std::vector<std::thread> threads;
    for (const auto& region : regions) {

        // Create a thread for each region, calling detectSVsFromRegion
        std::cout << "Detecting SVs in " << region << std::endl;
        std::cout << "Creating thread" << std::endl;
        std::thread t(&SVCaller::detectSVsFromRegion, this, region, std::ref(sv_calls), fp_in, bamHdr, idx, std::ref(this->mtx_sv_calls), std::ref(this->mtx_bam));
        threads.push_back(std::move(t));
    }

    // Join the threads
    for (auto& t : threads) {
        t.join();
    }
    std::cout << "Joined threads" << std::endl;

    // Close the BAM file
    sam_close(fp_in);

    // Destroy the header
    bam_hdr_destroy(bamHdr);

    // Destroy the index
    hts_idx_destroy(idx);

    // Print the number of SV calls
    std::cout << sv_calls.size() << " SV calls" << std::endl;

    // Return the SV calls
    return sv_calls;
}
