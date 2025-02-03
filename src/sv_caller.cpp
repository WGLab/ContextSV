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
#include <bitset>
#include <unordered_set>

#include "ThreadPool.h"
#include "utils.h"
#include "sv_types.h"
#include "version.h"
#include "fasta_query.h"
/// @endcond

# define DUP_SEQSIM_THRESHOLD 0.9  // Sequence similarity threshold for duplication detection

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
    // std::unordered_map<std::string, uint8_t> primary_map_qual;
    uint32_t num_alignments = 0;
    while (readNextAlignment(fp_in, itr, bam1) >= 0) {

        // Skip secondary and unmapped alignments, duplicates, QC failures, and low mapping quality
        if (bam1->core.flag & BAM_FSECONDARY || bam1->core.flag & BAM_FUNMAP || bam1->core.flag & BAM_FDUP || bam1->core.flag & BAM_FQCFAIL || bam1->core.qual < this->min_mapq) {
            continue;
        }
        const std::string qname = bam_get_qname(bam1);  // Query template name
        uint8_t mapq = bam1->core.qual;  // Mapping quality

        // Process primary alignments
        if (!(bam1->core.flag & BAM_FSUPPLEMENTARY)) {
            // primary_map[qname] = itr;
            // Store chromosome (TID), start, and end positions (1-based) of the
            // primary alignment, and the strand (true for forward, false for reverse)
            primary_map[qname] = GenomicRegion{bam1->core.tid, bam1->core.pos + 1, bam_endpos(bam1), !(bam1->core.flag & BAM_FREVERSE), mapq};
            // primary_map_qual[qname] = bam1->core.qual;
            primary_count++;

        // Process supplementary alignments
        } else if (bam1->core.flag & BAM_FSUPPLEMENTARY) {
            // supp_map[qname].push_back(itr);
            // Store chromosome (TID), start, and end positions (1-based) of the
            // supplementary alignment, and the strand (true for forward, false for reverse)
            supp_map[qname].push_back(GenomicRegion{bam1->core.tid, bam1->core.pos + 1, bam_endpos(bam1), !(bam1->core.flag & BAM_FREVERSE), mapq});
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
    printMessage(region + ": Found " + std::to_string(primary_map.size()) + " primary and " + std::to_string(supplementary_count) + " supplementary alignments");

    // Filter overlapping primary alignments and keep the one with the highest mapping
    // quality
    // std::vector<std::string> to_remove_overlapping;
    std::unordered_set<std::string> to_remove_overlapping;
    for (const auto& entry1 : primary_map) {
        const std::string& qname1 = entry1.first;
        const GenomicRegion& primary1 = entry1.second;
        for (const auto& entry2 : primary_map) {
            const std::string& qname2 = entry2.first;
            if (qname1 == qname2) {
                continue;
            }
            const GenomicRegion& primary2 = entry2.second;
            if (primary1.tid == primary2.tid && primary1.start <= primary2.end && primary1.end >= primary2.start) {
                // Overlapping primary alignments
                // printMessage("Overlapping primary alignments with quality " + std::to_string(primary_map_qual[qname1]) + " and " + std::to_string(primary_map_qual[qname2]));
                // if (primary_map_qual[qname1] < primary_map_qual[qname2]) {
                if (primary1.qual < primary2.qual) {
                    // to_remove_overlapping.push_back(qname1);
                    to_remove_overlapping.insert(qname1);
                } else {
                    // If equal, remove the shorter alignment
                    if (primary1.end - primary1.start < primary2.end - primary2.start) {
                        // to_remove_overlapping.push_back(qname1);
                        to_remove_overlapping.insert(qname1);
                    } else {
                        // to_remove_overlapping.push_back(qname2);
                        to_remove_overlapping.insert(qname2);
                    }
                }
            }
        }
    }

    for (const std::string& qname : to_remove_overlapping) {
        primary_map.erase(qname);
        supp_map.erase(qname);
    }
    printMessage(region + ": Removed " + std::to_string(to_remove_overlapping.size()) + " overlapping primary alignments");
    printMessage(region + ": Found " + std::to_string(primary_map.size()) + " primary and " + std::to_string(supp_map.size()) + " supplementary alignments after filtering");
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
    uint8_t qual = alignment->core.qual;

    // Loop through the CIGAR string, process operations, detect SVs (primary
    // only), and calculate sequence identity for potential duplications (primary only)
    uint32_t ref_pos;
    uint32_t ref_end;
    double default_lh = 0.0;
    const std::string amb_bases = "RYKMSWBDHV";  // Ambiguous bases
    std::bitset<256> amb_bases_bitset;
    for (char base : amb_bases) {
        amb_bases_bitset.set(base);
    }
    for (int i = 0; i < cigar_len; i++) {
        int op_len = bam_cigar_oplen(cigar[i]);  // CIGAR operation length
        int op = bam_cigar_op(cigar[i]);  // CIGAR operation
        if (op_len >= 50) {
            
            // Process the CIGAR operation
            if (op == BAM_CINS && is_primary) {

                // Get the sequence of the insertion from the query
                std::string ins_seq_str(op_len, ' ');
                for (int j = 0; j < op_len; j++) {
                    // Replace ambiguous bases with N
                    char base = seq_nt16_str[bam_seqi(bam_get_seq(alignment), query_pos + j)];
                    if (amb_bases_bitset.test(base)) {
                        ins_seq_str[j] = 'N';
                    } else {
                        ins_seq_str[j] = base;
                    }
                }
                
                // Before the insertion
                if (pos >= (uint32_t)op_len-1)
                {
                    uint32_t bp1 = pos - (op_len - 1) + 1;
                    uint32_t bp2 = bp1 + op_len - 1; //pos + 1;

                    if (ref_genome.compare(chr, bp1, bp2, ins_seq_str, DUP_SEQSIM_THRESHOLD))
                    {
                        int read_depth = this->calculateReadDepth(pos_depth_map, bp1, bp2);
                        addSVCall(sv_calls, bp1, bp2, SVType::DUP, "<DUP>", "LSEQSIM", "./.", default_lh, read_depth, qual);
                        continue;
                    }
                }

                // After the insertion
                if (pos + op_len < ref_genome.getChromosomeLength(chr))
                {
                    uint32_t bp1 = pos + 1;
                    uint32_t bp2 = bp1 + op_len - 1;

                    if (ref_genome.compare(chr, bp1, bp2, ins_seq_str, DUP_SEQSIM_THRESHOLD))
                    {
                        int read_depth = this->calculateReadDepth(pos_depth_map, bp1, bp2);
                        addSVCall(sv_calls, bp1, bp2, SVType::DUP, "<DUP>", "RSEQSIM", "./.", default_lh, read_depth, qual);
                        continue;
                    }
                }

                // Add as an insertion
                // For read depth calculation, use the previous and current
                // positions (1-based)
                uint32_t ins_pos = pos + 1;
                uint32_t ins_end = ins_pos + op_len - 1;
                int read_depth = this->calculateReadDepth(pos_depth_map, ins_pos-1, ins_pos);
                
                // Determine the ALT allele format based on small vs. large insertion
                std::string alt_allele = "<INS>";
                if (op_len <= 50) {
                    alt_allele = ins_seq_str;
                }
                
                addSVCall(sv_calls, ins_pos, ins_end, SVType::INS, alt_allele, "CIGARINS", "./.", default_lh, read_depth, qual);

            // Check if the CIGAR operation is a deletion
            } else if (op == BAM_CDEL && is_primary) {

                ref_pos = pos+1;
                ref_end = ref_pos + op_len -1;
                int read_depth = this->calculateReadDepth(pos_depth_map, ref_pos, ref_end);
                addSVCall(sv_calls, ref_pos, ref_end, SVType::DEL, "<DEL>", "CIGARDEL", "./.", default_lh, read_depth, qual);

                // Print if the ref pos is within the range 44007800-44007930
                if (ref_pos >= 44007800 && ref_pos <= 44007930) {
                    printMessage("DEL: " + chr + ":" + std::to_string(ref_pos) + "-" + std::to_string(ref_end) + " (LENGTH " + std::to_string(op_len) + ")");
                }
            }
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
    // int filter_threshold = 4;  // Minimum number of supporting reads for an SV call
    // int filter_threshold = 10;  // Minimum number of supporting reads for an
    // SV call
    int cigar_sv_support_threshold = input_data.getMinReadSupport();  // Minimum number of supporting reads for an SV call
    // int split_sv_support_threshold = 4;  // Minimum number of supporting
    // reads for an SV call
    int split_sv_support_threshold = input_data.getMinReadSupport();
    // printMessage("Processing chromosome " + chr + " with filter threshold: "
    // + std::to_string(filter_threshold));
    double dbscan_epsilon = input_data.getDBSCAN_Epsilon();
    int dbscan_min_pts = input_data.getDBSCAN_MinPts();

    // Open the BAM file
    std::string bam_filepath = input_data.getLongReadBam();
    samFile *fp_in = sam_open(bam_filepath.c_str(), "r");
    if (!fp_in) {
        printError("ERROR: failed to open " + bam_filepath);
        return;
    }

    // Set multi-threading
    int num_threads = input_data.getThreadCount();
    hts_set_threads(fp_in, num_threads);

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
    double mean_chr_cov = cnv_caller.calculateMeanChromosomeCoverage(chr, chr_pos_depth_map, bam_filepath, thread_count);
    if (mean_chr_cov == 0.0 || chr_pos_depth_map.size() == 0) {
        return;
    }

    // Detect SVs from the CIGAR strings
    printMessage(chr + ": CIGAR SVs...");
    this->detectCIGARSVs(fp_in, idx, bamHdr, region, chr_sv_calls, chr_pos_depth_map, ref_genome);

    printMessage(chr + ": Merging CIGAR...");
    // filterSVsWithLowSupport(chr_sv_calls, cigar_sv_support_threshold);
    mergeSVs(chr_sv_calls, dbscan_epsilon, dbscan_min_pts);
    // filterSVsWithLowSupport(chr_sv_calls, cigar_sv_support_threshold);
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
    this->detectSVsFromSplitReads(region, fp_in, idx, bamHdr, chr_sv_calls, cnv_caller, hmm, mean_chr_cov, chr_pos_depth_map, input_data);

    // Sort the SV calls by start position
    std::sort(chr_sv_calls.begin(), chr_sv_calls.end(), [](const SVCall& a, const SVCall& b) {
        return a.start < b.start;
    });

    // Merge the SV calls from the current region
    // printMessage(chr + ": Merging split reads...");
    // filterSVsWithLowSupport(chr_sv_calls, split_sv_support_threshold);
    // filterSVsWithLowSupport(chr_sv_calls, split_sv_support_threshold, "SPLIT");

    // Run a final merge on the combined SV calls
    // printMessage(chr + ": Merging final calls...");
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

    // Lambda to process a chromosome
    auto process_chr = [&](const std::string& chr) {
        try {
            std::vector<SVCall> sv_calls;
            InputData chr_input_data = input_data;  // Use a thread-local copy
            this->processChromosome(chr, hmm, sv_calls, chr_input_data, ref_genome);
            {
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
    int total_chr_count = futures.size();
    int current_chr = 0;
    for (auto& future : futures) {
        try {
            current_chr++;
            future.get();
            printMessage("Chromosome task "+ std::to_string(current_chr) + " of " + std::to_string(total_chr_count) + " completed.");
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
void SVCaller::detectSVsFromSplitReads(const std::string& region, samFile* fp_in, hts_idx_t* idx, bam_hdr_t* bamHdr, std::vector<SVCall>& sv_calls, const CNVCaller& cnv_caller, const CHMM& hmm, double mean_chr_cov, const std::vector<uint32_t>& pos_depth_map, const InputData& input_data)
{
    // printMessage(region + ": Getting split alignments...");
    std::unordered_map<std::string, GenomicRegion> primary_map;
    std::unordered_map<std::string, std::vector<GenomicRegion>> supp_map;
    this->getSplitAlignments(fp_in, idx, bamHdr, region, primary_map, supp_map);

    // Find split-read SV evidence
    // printMessage(region + ": Finding split-read SVs...");
    std::vector<SVCall> split_sv_calls;
    int current_primary = 0;
    int primary_count = primary_map.size();
    uint32_t min_cnv_length = input_data.getMinCNVLength();
    for (auto& entry : primary_map) {
        current_primary++;
        const std::string& qname = entry.first;
        GenomicRegion& primary = entry.second;
        const std::string& primary_chr = bamHdr->target_name[primary.tid];

      	// Find the largest supplementary alignment
        auto& supp_regions = supp_map[qname];
        // GenomicRegion largest_supp = supp_regions[0];
        auto it = std::max_element(supp_regions.begin(), supp_regions.end(), [](const GenomicRegion& a, const GenomicRegion& b) {
            return a.end - a.start < b.end - b.start;
        });
        GenomicRegion largest_supp = *it;

        // If on a different chromosome, label as a translocation
        if (primary.tid != largest_supp.tid) {
            // Note that these do not currently have a likelihood score or read depth
            // Create two BND records for the translocation
            // Create the alternate allele format for the first BND record
            const std::string& supp_chr = bamHdr->target_name[largest_supp.tid];
            std::string alt_allele = "N[" + supp_chr + ":" + std::to_string(largest_supp.start) + "[";
            if (largest_supp.strand == false) {
                // Reverse-oriented relative to the reference
                alt_allele = "N]" + supp_chr + ":" + std::to_string(largest_supp.start) + "]";
            }
            addSVCall(split_sv_calls, primary.start, primary.end, SVType::BND, alt_allele, "SPLIT", "./.", 0.0, 0, primary.qual);

            // Create the alternate allele format for the second BND record
            alt_allele = "N[" + primary_chr + ":" + std::to_string(primary.start) + "[";
            if (primary.strand == false) {
                // Reverse-oriented relative to the reference
                alt_allele = "N]" + primary_chr + ":" + std::to_string(primary.start) + "]";
            }
            addSVCall(split_sv_calls, largest_supp.start, largest_supp.end, SVType::BND, alt_allele, "SPLIT", "./.", 0.0, 0, largest_supp.qual);

            continue;
        }

        // Inversion detection
        bool is_opposite_strand = primary.strand != largest_supp.strand;
        if (is_opposite_strand) {
            // if (supp_length >= min_cnv_length) {
            if (largest_supp.end - largest_supp.start >= min_cnv_length) {

                // Print error if the start position is greater than the end
                // position
                // if (supp_start > supp_end) {
                if (largest_supp.start > largest_supp.end) {
                    printError("ERROR: Invalid inversion coordinates: " + primary_chr + ":" + std::to_string(largest_supp.start) + "-" + std::to_string(largest_supp.end));
                    // printError("ERROR: Invalid inversion coordinates: " + primary_chr + ":" + std::to_string(supp_start) + "-" + std::to_string(supp_end));
                    continue;
                }

                std::tuple<double, SVType, std::string, bool> result = cnv_caller.runCopyNumberPrediction(primary_chr, hmm, largest_supp.start, largest_supp.end, mean_chr_cov, pos_depth_map, input_data);
                if (std::get<1>(result) == SVType::UNKNOWN) {
                    continue;
                }

                double supp_lh = std::get<0>(result);
                SVType supp_type = std::get<1>(result);
                // printMessage("Test3");
                int read_depth = this->calculateReadDepth(pos_depth_map, largest_supp.start, largest_supp.end);
                // int read_depth = this->calculateReadDepth(pos_depth_map, supp_start, supp_end);
                if (supp_type == SVType::NEUTRAL) {
                    // addSVCall(sv_calls, supp_start, supp_end, "INV",
                    // "<INV>", "SPLIT", "./.", supp_lh, read_depth);
                    addSVCall(split_sv_calls, largest_supp.start, largest_supp.end, SVType::INV, "<INV>", "SPLIT", "./.", supp_lh, read_depth, largest_supp.qual);
                    continue;
                    
                } else if (supp_type == SVType::DUP) {
                    // addSVCall(sv_calls, supp_start, supp_end, "INVDUP",
                    // "<INV>", "SPLIT", "./.", supp_lh, read_depth);
                    addSVCall(split_sv_calls, largest_supp.start, largest_supp.end, SVType::INV_DUP, "<INV>", "SPLIT", "./.", supp_lh, read_depth, largest_supp.qual);
                    continue;
                }
            }
        }

        // Analyze split-read evidence for deletions and duplications
        uint8_t mean_qual = (primary.qual + largest_supp.qual) / 2;
        bool gap_exists = false;
        uint32_t boundary_left, boundary_right, gap_left, gap_right;
        boundary_left = std::min(primary.start, largest_supp.start);
        boundary_right = std::max(primary.end, largest_supp.end);
        gap_left = std::min(primary.end, largest_supp.start);
        gap_right = std::max(primary.start, largest_supp.end);
        gap_exists = gap_left < gap_right;
        
        // Run copy number variant predictions on the boundary if large enough
        if (boundary_right - boundary_left >= min_cnv_length) {

            // Print error if the start position is greater than the end
            // position
            if (boundary_left > boundary_right) {
                printError("ERROR: Invalid boundary coordinates: " + primary_chr + ":" + std::to_string(boundary_left) + "-" + std::to_string(boundary_right));
                continue;
            }

            // printMessage(region + ": Running copy number prediction for
            // boundary...");
            // printMessage("Running copy number prediction, length: " + std::to_string(boundary_right - boundary_left));
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

                // printMessage(region + ": Running copy number prediction for
                // gap...");
                // printMessage("Running copy number prediction, length: " + std::to_string(gap_right - gap_left));
                std::tuple<double, SVType, std::string, bool> gap_result = cnv_caller.runCopyNumberPrediction(primary_chr, hmm, gap_left, gap_right, mean_chr_cov, pos_depth_map, input_data);
                if (std::get<1>(gap_result) == SVType::UNKNOWN) {
                    continue;
                }
                double gap_lh = std::get<0>(gap_result);
                SVType gap_type = std::get<1>(gap_result);

                // If higher likelihood than the boundary, add the gap as the SV call
                if (gap_lh > bd_lh) {
                    int read_depth = this->calculateReadDepth(pos_depth_map, gap_left, gap_right);
                    std::string alt_allele = gap_type == SVType::NEUTRAL ? "." : "<" + getSVTypeString(gap_type) + ">";
                    addSVCall(split_sv_calls, gap_left, gap_right, gap_type, alt_allele, "SPLIT", "./.", gap_lh, read_depth, mean_qual);
                } else {
                    // Add the boundary as the SV call
                    int read_depth = this->calculateReadDepth(pos_depth_map, boundary_left, boundary_right);
                    std::string alt_allele = bd_type == SVType::NEUTRAL ? "." : "<" + getSVTypeString(bd_type) + ">";
                    addSVCall(split_sv_calls, boundary_left, boundary_right, bd_type, alt_allele, "SPLIT", "./.", bd_lh, read_depth, mean_qual);
                }
            } else {
                // Add the boundary as the SV call
                int read_depth = this->calculateReadDepth(pos_depth_map, boundary_left, boundary_right);
                std::string alt_allele = bd_type == SVType::NEUTRAL ? "." : "<" + getSVTypeString(bd_type) + ">";
                addSVCall(split_sv_calls, boundary_left, boundary_right, bd_type, alt_allele, "SPLIT", "./.", bd_lh, read_depth, mean_qual);
            }
        }

        // Print progress every 1000 primary alignments
        if (current_primary % 1000 == 0) {
            printMessage(region + ": Processed " + std::to_string(current_primary) + " of " + std::to_string(primary_count) + " primary alignments...");
        }
    }

    // Merge the split-read SV calls
    printMessage(region + ": Merging split-read SVs...");
    mergeSVs(split_sv_calls, 0.1, 2);

    // Unify the SV calls
    sv_calls.insert(sv_calls.end(), split_sv_calls.begin(), split_sv_calls.end());
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
        "##INFO=<ID=CLUSTER,Number=1,Type=Integer,Description=\"Cluster size\">",
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
            std::string sv_type_str = getSVTypeString(sv_call.sv_type);
            std::string genotype = sv_call.genotype;
            std::string data_type_str = sv_call.data_type;
            std::string alt_allele = sv_call.alt_allele;
            double hmm_likelihood = sv_call.hmm_likelihood;
            int sv_length = end - start + 1;
            int cluster_size = sv_call.cluster_size;
            /*
            if (sv_type_str == "DEL") {
            	sv_length++;
        	}
        	*/
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

            // Print the REF allele if SVTYPE = DUP and if it is empty or "." (symbolic)
            if (sv_type_str == "DUP" && (ref_allele == "" || ref_allele == ".")) {
                printMessage("REF allele for DUP at " + chr + ":" + std::to_string(start) + "-" + std::to_string(end) + ": " + ref_allele + ", ALT allele: " + alt_allele);
            }

            // Create the VCF parameter strings
            std::string info_str = "END=" + std::to_string(end) + ";SVTYPE=" + sv_type_str + \
                ";SVLEN=" + std::to_string(sv_length) + ";SVMETHOD=" + sv_method + ";ALN=" + data_type_str + \
                ";HMM=" + std::to_string(hmm_likelihood) + ";SUPPORT=" + std::to_string(support) + ";CLUSTER=" + std::to_string(cluster_size);
                
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
    vcf_stream.close();

    std::cout << "Saved SV calls to " << output_vcf << std::endl;

    // Create a compressed and indexed VCF file
    std::cout << "Creating compressed and indexed VCF file..." << std::endl;
    std::string bgzip_cmd = "bgzip -f " + output_vcf;
    std::string tabix_cmd = "tabix -p vcf " + output_vcf + ".gz";
    std::system(bgzip_cmd.c_str());
    std::system(tabix_cmd.c_str());

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
