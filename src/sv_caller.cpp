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

#include "utils.h"
/// @endcond

int SVCaller::readNextAlignment(samFile *fp_in, hts_itr_t *itr, bam1_t *bam1)
{
    // Lock the mutex while reading the next alignment
//    std::lock_guard<std::mutex> lock(this->bam_mtx);

    // Read the next alignment
    int ret = sam_itr_next(fp_in, itr, bam1);

    // Return the result of reading the next alignment
    return ret;
}

// SVData SVCaller::detectSVsFromRegion(std::string region, samFile *fp_in, bam_hdr_t *bamHdr, hts_idx_t *idx)
SVData SVCaller::detectSVsFromRegion(std::string region)
{
    SVData sv_calls;
    std::string bam_filepath = this->input_data->getLongReadBam();

    // Lock the mutex while setting up the BAM file
//    this->bam_mtx.lock();

    // Check if the CIGAR string should be used for SV calling
    bool disable_cigar = this->input_data->getDisableCIGAR();

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

    // Unlock the mutex after setting up the BAM file
//    this->bam_mtx.unlock();

    // Loop through the alignments
    // Create a map of primary and supplementary alignments by QNAME (query template name)
    int num_alignments = 0;
    QueryMap primary_alignments;  // TODO: Add depth to primary alignments
    QueryMap supplementary_alignments;

    //std::cout << "Reading alignments..." << std::endl;
    // Thread-safe printing
    printMessage("Detecting SVs from " + region);

    //printMessage("Starting SV count: " + std::to_string(sv_calls.size()));
    while (readNextAlignment(fp_in, itr, bam1) >= 0) {

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
//                this->bam_mtx.lock();
                std::string chr = bamHdr->target_name[bam1->core.tid];
                int64_t start = bam1->core.pos;
                int64_t end = bam_endpos(bam1);
//                this->bam_mtx.unlock();

                AlignmentData alignment(chr, start, end, ".");
                primary_alignments[qname].push_back(alignment);
            
                // Call SVs directly from the CIGAR string
                if (disable_cigar == false) {

                    // printMessage("Calling CIGAR SVs from " + region, print_mtx);
                    this->detectSVsFromCIGAR(bamHdr, bam1, sv_calls);
                    // printMessage("Finished calling CIGAR SVs from " + region, print_mtx);
                }

            // Process supplementary alignments
            } else if (bam1->core.flag & BAM_FSUPPLEMENTARY) {

                // Add the supplementary alignment to the map
//                this->bam_mtx.lock();
                std::string chr = bamHdr->target_name[bam1->core.tid];
                int32_t start = bam1->core.pos;
                int32_t end = bam_endpos(bam1);
                AlignmentData alignment(chr, start, end, ".");
//                this->bam_mtx.unlock();

                // Add the supplementary alignment to the map
                supplementary_alignments[qname].push_back(alignment);
            }
        }

        // Increment the number of alignment records processed
        num_alignments++;

        // Print the number of alignments processed every 10 thousand
        if (num_alignments % 10000 == 0) {
            std::cout << "Region: " << region << " Alignments processed: " << num_alignments << std::endl;
        }
    }
    
    // Print the number of alignments processed
    //std::cout << num_alignments << " alignments processed" << std::endl;

    // Loop through the map of primary alignments by QNAME and find gaps and
    // overlaps from supplementary alignments
    //std::cout << "Running split read analysis..." << std::endl;
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
                //std::cout << "Supplementary alignment on different chromosome" << std::endl;
                continue;
            }

            // Get the start and end positions of the supplementary alignment
            int32_t supp_start = std::get<1>(supp_alignment);
            int32_t supp_end = std::get<2>(supp_alignment);

            // Gap analysis (deletion or duplication)
            if (supp_start < primary_start && supp_end < primary_start) {
                // Gap with supplementary before primary:
                // [supp_start] [supp_end] -- [primary_start] [primary_end]

                // Use the gap ends as the SV endpoints
                if (primary_start - supp_end >= this->min_sv_size) {
                    sv_calls.add(supp_chr, supp_end+1, primary_start+1, SVData::UNKNOWN, ".", "GAPINNER_A");
                }

                // ALso use the alignment ends as the SV endpoints
                if (primary_end - supp_start >= this->min_sv_size) {
                    sv_calls.add(supp_chr, supp_start+1, primary_end+1, SVData::UNKNOWN, ".", "GAPOUTER_A");
                }

                
            } else if (supp_start > primary_end && supp_end > primary_end) {
                // Gap with supplementary after primary:
                // [primary_start] [primary_end] -- [supp_start] [supp_end]

                // Use the gap ends as the SV endpoints
                if (supp_start - primary_end >= this->min_sv_size) {
                    sv_calls.add(supp_chr, primary_end+1, supp_start+1, SVData::UNKNOWN, ".", "GAPINNER_B");
                }

                // Also use the alignment ends as the SV endpoints
                if (supp_end - primary_start >= this->min_sv_size) {
                    sv_calls.add(supp_chr, primary_start+1, supp_end+1, SVData::UNKNOWN, ".", "GAPOUTER_B");
                }
            }
        }
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

    if (sv_calls.size() > 0) {
        printMessage("Found " + std::to_string(sv_calls.size()) + " SVs from " + region);
    }
    // printMessage("Finished detecting " + std::to_string(sv_calls.size()) + " SVs from " + region);

    // Return the SV calls
    return sv_calls;

    // Print the number of SV calls
    //std::cout << sv_calls.size() << " SV calls from " << region << std::endl;
}

SVCaller::SVCaller(InputData &input_data)
{
    this->input_data = &input_data;
}

// Main function for SV detection
SVData SVCaller::run()
{
    // Get SV calls from split read alignments (primary and supplementary) and
    // directly from the CIGAR string
    SVData sv_calls = this->detectSVsFromSplitReads();

    return sv_calls;
}

// Detect SVs from the CIGAR string of a read alignment.
void SVCaller::detectSVsFromCIGAR(bam_hdr_t* header, bam1_t* alignment, SVData& sv_calls)
{
//    this->bam_mtx.lock();

    // Get the chromosome
    std::string chr = header->target_name[alignment->core.tid];

    // Get the position
    int32_t pos = alignment->core.pos;

    // Get the CIGAR string
    uint32_t* cigar = bam_get_cigar(alignment);

    // Get the CIGAR length
    int cigar_len = alignment->core.n_cigar;

//    this->bam_mtx.unlock();

    // Track the query position
    int query_pos = 0;

    // Loop through the CIGAR string (0-based) and detect insertions and deletions in
    // reference coordinates (1-based)
    // POS is the leftmost position of where the alignment maps to the reference:
    // https://genome.sph.umich.edu/wiki/SAM

    // Get the reference genome
    // FASTAQuery ref_genome = this->input_data->getRefGenome();

    // This will be done using threading. The CIGAR string will be split into
    // chunks and each chunk will be processed in a separate thread. The
    // resulting SV calls will be concatenated at the end.

    std::vector<std::thread> threads;
    std::vector<SVData> sv_calls_vec;

    // Loop through the CIGAR string and process operations
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

                // To determine whether the insertion is a duplication, check
                // for sequence identity between the insertion and the
                // reference genome (> 90% identity is likely a duplication).

                // Loop starting from the reference position minus the
                // insertion length to the reference position plus the
                // insertion length, and obtain the highest sequence identity
                // using a sliding window of the insertion length
                bool is_duplication = false;

                // Tricky region START
                // Get the reference sequence for the region (+/- insertion length)
//                this->bam_mtx.lock();
                // this->query_mtx.lock();
                // std::string ref_seq_str = this->input_data->getRefGenome().query(chr, pos - op_len, pos + op_len - 1);
                // this->query_mtx.unlock();

                // // Check if the insertion sequence before or after the insertion position
                // // matches the reference sequence with > 50% identity
                // float target_seq_identity = 0.5;
                // int ref_seq_start = 0;
                // int ref_seq_end = (int) ref_seq_str.length() - (int) op_len;
                // std::string ref_seq_before = ref_seq_str.substr(ref_seq_start, op_len);
                // std::string ref_seq_after = ref_seq_str.substr(ref_seq_end, op_len);
                // int num_matches_before = 0;
                // int num_matches_after = 0;
                // for (int j = 0; j < op_len; j++) {
                //     if (ins_seq_str[j] == ref_seq_before[j]) {
                //         num_matches_before++;
                //     }
                //     if (ins_seq_str[j] == ref_seq_after[j]) {
                //         num_matches_after++;
                //     }
                // }
                // float seq_identity_before = (float)num_matches_before / (float)op_len;
                // float seq_identity_after = (float)num_matches_after / (float)op_len;
                // if (seq_identity_before > target_seq_identity || seq_identity_after > target_seq_identity) {
                //     is_duplication = true;
                // }
                // Tricky region END

                // Add to SV calls (1-based) with the appropriate SV type
                ref_pos = pos+1;
                ref_end = ref_pos + op_len -1;
                if (is_duplication) {
                    sv_calls.add(chr, ref_pos, ref_end, SVData::DUP, ins_seq_str, "CIGARDUP");
                    //std::cout << "ADDED CIGAR SV" << std::endl;
                } else {
                    sv_calls.add(chr, ref_pos, ref_end, SVData::INS, ins_seq_str, "CIGARINS");
                    //std::cout << "ADDED CIGAR SV" << std::endl;
                }
            }

        // Check if the CIGAR operation is a deletion
        } else if (op == BAM_CDEL) {

            // Add the SV if greater than the minimum SV size
            if (op_len >= this->min_sv_size) {
                
                // Add the deletion to the SV calls (1-based)
                ref_pos = pos+1;
                ref_end = ref_pos + op_len -1;
                sv_calls.add(chr, ref_pos, ref_end, SVData::DEL, ".", "CIGARDEL");
                //std::cout << "ADDED CIGAR SV" << std::endl;
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
SVData SVCaller::detectSVsFromSplitReads()
{
    // Open the BAM file
    std::string bam_filepath = this->input_data->getLongReadBam();

    // Get the region data
    bool whole_genome = this->input_data->getWholeGenome();
    std::vector<std::string> regions;
    if (whole_genome) {
        regions = this->input_data->getRefGenome().getChromosomes();
    } else {
        regions.push_back(this->input_data->getRegion());
    }

    // ----- Non-threading -----
    // Loop through each region and detect SVs
    std::cout << "Detecting SVs from " << regions.size() << " regions..." << std::endl;
    auto start1 = std::chrono::high_resolution_clock::now();
    // std::vector<SVData> sv_calls_vec;
    SVData sv_calls;
    for (const auto& region : regions) {
        // SVData sv_calls = this->detectSVsFromRegion(region);
        // sv_calls_vec.push_back(sv_calls);

        // Split the region into chunks and process each chunk in a separate
        // thread
        // Get the number of available threads
        int num_threads = std::thread::hardware_concurrency();

        // Get the chromosome length
        int chr_len = this->input_data->getRefGenome().getChromosomeLength(region);

        // Split the region into chunks
        int chunk_size = chr_len / num_threads;
        std::vector<std::string> region_chunks;
        for (int i = 0; i < num_threads; i++) {
            int start = i * chunk_size;
            int end = start + chunk_size;
            if (i == num_threads - 1) {
                end = chr_len;
            }
            std::string chunk = region + ":" + std::to_string(start) + "-" + std::to_string(end);
            region_chunks.push_back(chunk);
        }

        // Loop through the chunks and process each chunk in a separate thread
        std::vector<SVData> sv_calls_vec;
        std::vector<std::thread> threads;

        // Vector of futures
        std::vector<std::future<SVData>> futures;
        for (const auto& sub_region : region_chunks) {
            // std::cout << "Sub-region: " << sub_region << std::endl;

            // Detect SVs from the sub-region in a separate thread using a
            // future
            std::future<SVData> future = std::async(std::launch::async, &SVCaller::detectSVsFromRegion, this, sub_region);

            // Add the future to the list of futures
            futures.push_back(std::move(future));
        }

        // Get the SV calls from the futures
        for (auto& future : futures) {
            // Wait for the future to be ready
            future.wait();

            // Get the SV calls from the future
            SVData sv_calls_region = future.get();
            //std::cout << "Future SV calls: " << sv_calls_region.size() << std::endl;
            sv_calls_vec.push_back(sv_calls_region);
        }

        // Combine the SV calls from each region
        std::cout << "Combining SV calls..." << std::endl;
        auto start2 = std::chrono::high_resolution_clock::now();
        for (auto it = sv_calls_vec.begin(); it != sv_calls_vec.end(); ++it) {
            //std::cout << "SV calls: " << it->size() << std::endl;
            sv_calls.concatenate(*it);
            //std::cout << "Updated SV calls: " << sv_calls.size() << std::endl;
        }

        auto end2 = std::chrono::high_resolution_clock::now();
        std::cout << "Combining SV calls finished. Elapsed time: " << getElapsedTime(start2, end2) << std::endl;
    }
    auto end1 = std::chrono::high_resolution_clock::now();
    std::cout << "Finished detecting SVs from " << regions.size() << " regions. Elapsed time (h:m:s): " << getElapsedTime(start1, end1) << std::endl;

    // Return the SV calls
    return sv_calls;
}
