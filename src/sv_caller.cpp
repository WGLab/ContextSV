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
#include "dbscan.h"
#include "dbscan1d.h"
/// @endcond

# define DUP_SEQSIM_THRESHOLD 0.9  // Sequence similarity threshold for duplication detection

//std::mutex bam_mutex;

int SVCaller::readNextAlignment(samFile *fp_in, hts_itr_t *itr, bam1_t *bam1)
{
    std::lock_guard<std::mutex> lock(this->shared_mutex);
    int ret = sam_itr_next(fp_in, itr, bam1);
    return ret;
}

std::vector<std::string> SVCaller::getChromosomes(const std::string &bam_filepath)
{
    // Open the BAM file
    samFile *fp_in = sam_open(bam_filepath.c_str(), "r");
    if (!fp_in) {
        printError("ERROR: failed to open BAM file " + bam_filepath);
        return {};
    }
    bam_hdr_t *bamHdr = sam_hdr_read(fp_in);
    if (!bamHdr) {
        sam_close(fp_in);
        printError("ERROR: failed to read header from " + bam_filepath);
        return {};
    }
    std::vector<std::string> chromosomes;
    for (int i = 0; i < bamHdr->n_targets; i++) {
        chromosomes.push_back(bamHdr->target_name[i]);
        // printMessage("Chromosome: " + std::string(bamHdr->target_name[i]));
    }
    bam_hdr_destroy(bamHdr);
    sam_close(fp_in);
    return chromosomes;
}

std::vector<SVCall> SVCaller::getSplitAlignments(samFile *fp_in, hts_idx_t *idx, bam_hdr_t *bamHdr, const std::string &region, std::unordered_map<std::string, GenomicRegion> &primary_map, std::unordered_map<std::string, std::vector<GenomicRegion>> &supp_map)
{
    // Create a read and iterator for the region
    bam1_t *bam1 = bam_init1();
    if (!bam1) {
        printError("ERROR: failed to initialize BAM record");
        return {};
    }
    hts_itr_t *itr = sam_itr_querys(idx, bamHdr, region.c_str());
    if (!itr) {
        bam_destroy1(bam1);
        printError("ERROR: failed to query region " + region);
        return {};
    }

    uint32_t primary_count = 0;
    uint32_t supplementary_count = 0;

    // Main loop to process the alignments
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
            // Store chromosome (TID), start, and end positions (1-based) of the
            // primary alignment, and the strand (true for forward, false for reverse)
            primary_map[qname] = GenomicRegion{bam1->core.tid, bam1->core.pos + 1, bam_endpos(bam1), !(bam1->core.flag & BAM_FREVERSE), mapq, 0};
            primary_count++;

        // Process supplementary alignments
        } else if (bam1->core.flag & BAM_FSUPPLEMENTARY) {
            // Store chromosome (TID), start, and end positions (1-based) of the
            // supplementary alignment, and the strand (true for forward, false for reverse)
            supp_map[qname].push_back(GenomicRegion{bam1->core.tid, bam1->core.pos + 1, bam_endpos(bam1), !(bam1->core.flag & BAM_FREVERSE), mapq, 0});
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

    // Identify overlapping primary alignments and then cluster their primary
    // start, end vs. supplementary alignment start, end positions, keeping the
    // median of the largest cluster for the primary and supplementary positions
    // as the final genome coordinates of the SV
    IntervalNode* root = nullptr;
    for (const auto& entry : primary_map) {
        const std::string& qname = entry.first;
        const GenomicRegion& region = entry.second;
        root = insert(root, region, qname);
    }
    std::vector<std::vector<std::string>> primary_clusters;
    std::set<std::string> processed;

    for (const auto& entry : primary_map) {
        const std::string& qname = entry.first;
        if (processed.find(qname) != processed.end()) {
            continue;  // Skip already processed primary alignments
        }
        const GenomicRegion& region = entry.second;
        std::vector<std::string> overlap_group;
        findOverlaps(root, region, overlap_group);
        for (const std::string& qname : overlap_group) {
            processed.insert(qname);
        }
        if (overlap_group.size() > 1) {
            primary_clusters.push_back(overlap_group);
        }
    }
    printMessage(region + ": Found " + std::to_string(primary_clusters.size()) + " groups of overlapping primary alignments");

    // For each primary alignment cluster the supplementary alignment start and
    // end positions, keeping the median of the largest cluster
    std::vector<SVCall> sv_candidates;
    int current_group = 0;
    int min_length = 2000;
    int max_length = 1000000;
    for (const auto& primary_group : primary_clusters) {
        // Determine if the primary alignments are mostly on opposite strands to
        // the corresponding supplementary alignments (potential inversions)
        bool inversion = false;
        for (const std::string& qname : primary_group) {
            const std::vector<GenomicRegion>& regions = supp_map[qname];
            int num_supp = (int) regions.size();
            int num_opposite_strand = 0;
            for (const GenomicRegion& region : regions) {
                if (region.strand != primary_map[qname].strand) {
                    num_opposite_strand++;
                }
            }
            if (static_cast<double>(num_opposite_strand) / static_cast<double>(num_supp) > 0.5) {
                inversion = true;
            }
        }

        // Use DBSCAN to cluster primary alignment start, end positions
        DBSCAN1D dbscan(100, 5);
        current_group++;
        std::vector<int> starts;
        std::vector<int> ends;
        std::vector<bool> primary_strands;
        for (const std::string& qname : primary_group) {
            const GenomicRegion& region = primary_map[qname];
            starts.push_back(region.start);
            ends.push_back(region.end);
            primary_strands.push_back(region.strand);
        }

        // Get the largest cluster of primary alignment start positions
        dbscan.fit(starts);
        std::vector<int> primary_start_cluster = dbscan.getLargestCluster(starts);

        // Get the largest cluster of primary alignment end positions
        dbscan.fit(ends);
        std::vector<int> primary_end_cluster = dbscan.getLargestCluster(ends);

        // Continue if no clusters were found
        if (primary_start_cluster.empty() && primary_end_cluster.empty()) {
            continue;
        }

        // Get the supplementary alignment positions
        std::vector<int> supp_starts;
        std::vector<int> supp_ends;
        std::vector<bool> supp_strands;
        for (const std::string& qname : primary_group) {
            const std::vector<GenomicRegion>& regions = supp_map[qname];
            for (const GenomicRegion& region : regions) {
                supp_starts.push_back(region.start);
                supp_ends.push_back(region.end);
                supp_strands.push_back(region.strand);
            }
        }

        // Get the largest cluster of supplementary alignment start positions
        dbscan.fit(supp_starts);
        std::vector<int> supp_start_cluster = dbscan.getLargestCluster(supp_starts);

        // Get the largest cluster of supplementary alignment end positions
        dbscan.fit(supp_ends);
        std::vector<int> supp_end_cluster = dbscan.getLargestCluster(supp_ends);

        // Continue if no clusters were found
        if (supp_start_cluster.empty() && supp_end_cluster.empty()) {
            continue;
        }

        // Use the median of the largest cluster of primary and supplementary
        // alignment start, end positions as the final genome coordinates of the
        // SV
        int primary_pos = -1;
        int primary_pos2 = -1;
        if (primary_start_cluster.size() > primary_end_cluster.size()) {
            std::sort(primary_start_cluster.begin(), primary_start_cluster.end());
            primary_pos = primary_start_cluster[primary_start_cluster.size() / 2];
        } else if (primary_end_cluster.size() > primary_start_cluster.size()) {
            std::sort(primary_end_cluster.begin(), primary_end_cluster.end());
            primary_pos = primary_end_cluster[primary_end_cluster.size() / 2];
        } else {
            // Use both positions
            std::sort(primary_start_cluster.begin(), primary_start_cluster.end());
            std::sort(primary_end_cluster.begin(), primary_end_cluster.end());
            primary_pos = primary_start_cluster[primary_start_cluster.size() / 2];
            primary_pos2 = primary_end_cluster[primary_end_cluster.size() / 2];
        }

        // Get the supplementary alignment positions
        int supp_pos = -1;
        int supp_pos2 = -1;
        if (supp_start_cluster.size() > supp_end_cluster.size()) {
            std::sort(supp_start_cluster.begin(), supp_start_cluster.end());
            supp_pos = supp_start_cluster[supp_start_cluster.size() / 2];
        } else if (supp_end_cluster.size() > supp_start_cluster.size()) {
            std::sort(supp_end_cluster.begin(), supp_end_cluster.end());
            supp_pos = supp_end_cluster[supp_end_cluster.size() / 2];
        } else {
            // Use both positions. This has been shown to occur in nested SVs
            std::sort(supp_start_cluster.begin(), supp_start_cluster.end());
            std::sort(supp_end_cluster.begin(), supp_end_cluster.end());
            supp_pos = supp_start_cluster[supp_start_cluster.size() / 2];
            supp_pos2 = supp_end_cluster[supp_end_cluster.size() / 2];
        }

        // If two of either were found, use the larger SV candidate
        if (primary_pos2 != -1) {
            int sv_length1 = std::abs(primary_pos - supp_pos);
            int sv_length2 = std::abs(primary_pos2 - supp_pos);
            if (sv_length2 > sv_length1) {
                primary_pos = primary_pos2;
            }
        }
        if (supp_pos2 != -1) {
            int sv_length1 = std::abs(primary_pos - supp_pos);
            int sv_length2 = std::abs(primary_pos - supp_pos2);
            if (sv_length2 > sv_length1) {
                supp_pos = supp_pos2;
            }
        }

        if (primary_pos == -1 || supp_pos == -1) {
            continue;
        }

        // Store the SV candidate if the length is within the specified range
        int sv_start = std::min(primary_pos, supp_pos);
        int sv_end = std::max(primary_pos, supp_pos);
        int sv_length = sv_end - sv_start + 1;
        SVType sv_type = inversion ? SVType::INV : SVType::UNKNOWN;
        if (sv_length >= min_length && sv_length <= max_length) {
            SVCall sv_candidate(sv_start, sv_end, sv_type, ".", "NA", "./.", 0.0, 0, 0, 0);
            sv_candidates.push_back(sv_candidate);
        }
    }

    return sv_candidates;
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

    // Loop through the CIGAR string, process operations, detect SVs (primary
    // only), and calculate sequence identity for potential duplications (primary only)
    uint32_t ref_pos;
    uint32_t ref_end;
    double default_lh = 0.0;
    const std::string amb_bases = "RYKMSWBDHV";  // Ambiguous bases
    std::bitset<256> amb_bases_bitset;
    for (char base : amb_bases) {
        amb_bases_bitset.set(base);
        amb_bases_bitset.set(std::tolower(base));
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
                        // addSVCall(sv_calls, bp1, bp2, SVType::DUP, "<DUP>",
                        // "LSEQSIM", "./.", default_lh, read_depth);
                        SVCall sv_call(bp1, bp2, SVType::DUP, "<DUP>", "LSEQSIM", "./.", default_lh, read_depth, 1, 0);
                        addSVCall(sv_calls, sv_call);
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
                        SVCall sv_call(bp1, bp2, SVType::DUP, "<DUP>", "RSEQSIM", "./.", default_lh, read_depth, 1, 0);
                        addSVCall(sv_calls, sv_call);
                        // addSVCall(sv_calls, bp1, bp2, SVType::DUP, "<DUP>", "RSEQSIM", "./.", default_lh, read_depth);
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
                SVCall sv_call(ins_pos, ins_end, SVType::INS, alt_allele, "CIGARINS", "./.", default_lh, read_depth, 1, 0);
                addSVCall(sv_calls, sv_call);                
                // addSVCall(sv_calls, ins_pos, ins_end, SVType::INS, alt_allele, "CIGARINS", "./.", default_lh, read_depth);

            // Check if the CIGAR operation is a deletion
            } else if (op == BAM_CDEL && is_primary) {

                ref_pos = pos+1;
                ref_end = ref_pos + op_len -1;
                int read_depth = this->calculateReadDepth(pos_depth_map, ref_pos, ref_end);
                // addSVCall(sv_calls, ref_pos, ref_end, SVType::DEL, "<DEL>",
                // "CIGARDEL", "./.", default_lh, read_depth);
                SVCall sv_call(ref_pos, ref_end, SVType::DEL, "<DEL>", "CIGARDEL", "./.", default_lh, read_depth, 1, 0);
                addSVCall(sv_calls, sv_call);

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
    // uint32_t chr_len = bamHdr->target_len[bam_name2id(bamHdr, chr.c_str())];
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

    // Calculate the median read depth across the SV calls
    printMessage(chr + ": Calculating median SV read depth...");
    uint32_t cumulative_depth = 0;
    for (auto& sv_call : chr_sv_calls) {
        cumulative_depth += sv_call.read_depth;
    }
    double median_sv_depth = (double)cumulative_depth / (double)chr_sv_calls.size();
    printMessage("Median SV read depth: " + std::to_string(median_sv_depth));

    printMessage(chr + ": Merging CIGAR...");
    mergeSVs(chr_sv_calls, dbscan_epsilon, dbscan_min_pts);

    int region_sv_count = getSVCount(chr_sv_calls);
    printMessage("Total SVs detected from CIGAR string: " + std::to_string(region_sv_count));

    // Run copy number variant predictions on the SVs detected from the
    // CIGAR string, using a minimum CNV length threshold
    if (region_sv_count > 0) {
        printMessage(chr + ": CIGAR predictions...");
        cnv_caller.runCIGARCopyNumberPrediction(chr, chr_sv_calls, hmm, mean_chr_cov, chr_pos_depth_map, input_data);
    }

    // Run split-read SV and copy number variant predictions
    printMessage(chr + ": Split read SVs...");
    std::vector<SVCall> split_sv_calls;
    this->detectSVsFromSplitReads(region, fp_in, idx, bamHdr, split_sv_calls, cnv_caller, hmm, mean_chr_cov, chr_pos_depth_map, input_data);

    // // Merge the split-read SVs separately
    printMessage(chr + ": Merging split reads...");
    double split_epsilon = 0.45;
    int split_min_pts = 2;  // This is low since split alignments were already previously merged
    mergeSVs(split_sv_calls, split_epsilon, split_min_pts);

    printMessage(chr + ": Unifying SVs...");
    chr_sv_calls.insert(chr_sv_calls.end(), split_sv_calls.begin(), split_sv_calls.end());

    // mergeSVSubsets(chr_sv_calls);

    // Sort the SV calls by start position
    std::sort(chr_sv_calls.begin(), chr_sv_calls.end(), [](const SVCall& a, const SVCall& b) {
        return a.start < b.start;
    });

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
        // Get the chromosome from the user input argument
        chromosomes.push_back(input_data.getChromosome());
    } else {
        // chromosomes = ref_genome.getChromosomes();
        // Get the chromosomes from the input BAM file
        chromosomes = this->getChromosomes(input_data.getLongReadBam());
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
void SVCaller::detectSVsFromSplitReads(const std::string& region, samFile* fp_in, hts_idx_t* idx, bam_hdr_t* bamHdr, std::vector<SVCall>& split_sv_calls, const CNVCaller& cnv_caller, const CHMM& hmm, double mean_chr_cov, const std::vector<uint32_t>& pos_depth_map, const InputData& input_data)
{
    // printMessage(region + ": Getting split alignments...");
    std::unordered_map<std::string, GenomicRegion> primary_map;
    std::unordered_map<std::string, std::vector<GenomicRegion>> supp_map;
    std::vector<SVCall> sv_candidates = this->getSplitAlignments(fp_in, idx, bamHdr, region, primary_map, supp_map);

    // Run copy number predictions on the SVs detected from the split reads
    printMessage(region + ": Split read predictions...");
    int current_sv = 0;
    int total_svs = sv_candidates.size();
    for (auto& sv_candidate : sv_candidates) {
        bool is_inversion = sv_candidate.sv_type == SVType::INV;

        std::tuple<double, SVType, std::string, bool> result = cnv_caller.runCopyNumberPrediction(region, hmm, sv_candidate.start, sv_candidate.end, mean_chr_cov, pos_depth_map, input_data);
        if (std::get<1>(result) == SVType::UNKNOWN) {
            continue;
        }

        double supp_lh = std::get<0>(result);
        SVType supp_type = std::get<1>(result);
        std::string genotype = std::get<2>(result);

        // int read_depth = this->calculateReadDepth(pos_depth_map, sv_candidate.start, sv_candidate.end);
        if (supp_type != SVType::UNKNOWN) {
            if (is_inversion) {
                if (supp_type == SVType::DEL) {
                    supp_type = SVType::INV_DEL;
                } else if (supp_type == SVType::DUP) {
                    supp_type = SVType::INV_DUP;
                } else if (supp_type == SVType::NEUTRAL) {
                    supp_type = SVType::INV;
                }
            }
            
            if (supp_type != SVType::NEUTRAL) {
                int read_depth = this->calculateReadDepth(pos_depth_map, sv_candidate.start, sv_candidate.end);
                std::string alt_allele = "<" + getSVTypeString(supp_type) + ">";
                SVCall sv_call(sv_candidate.start, sv_candidate.end, supp_type, alt_allele, "SPLIT", genotype, supp_lh, read_depth, 1, sv_candidate.cluster_size);
                addSVCall(split_sv_calls, sv_call);
            }
        }
        current_sv++;
        if (current_sv % 1000 == 0) {
            printMessage("Processed " + std::to_string(current_sv) + " of " + std::to_string(total_svs) + " SV candidates");
        }
    }
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
        "##INFO=<ID=SVTYPE2,Number=1,Type=String,Description=\"Type of structural variant (if more than one)\">",
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
    int total_count = 0;
    for (const auto& pair : sv_calls) {
        std::string chr = pair.first;
        const std::vector<SVCall>& sv_calls = pair.second;
        std::cout << "Saving SV calls for " << chr << "..." << std::endl;
        for (const auto& sv_call : sv_calls) {
            // Get the SV candidate and SV info
            uint32_t start = sv_call.start;
            uint32_t end = sv_call.end;
            SVType sv_type = sv_call.sv_type;
            std::string genotype = sv_call.genotype;
            std::string data_type_str = sv_call.data_type;
            std::string alt_allele = sv_call.alt_allele;
            double hmm_likelihood = sv_call.hmm_likelihood;
            int sv_length = end - start + 1;
            int cluster_size = sv_call.cluster_size;
            int read_depth = sv_call.read_depth;
            std::string ref_allele = ".";
            int support = sv_call.support;

            // If the SV type is unknown, print a warning and skip
            if (sv_type == SVType::UNKNOWN || sv_type == SVType::NEUTRAL) {
                std::cerr << "Warning: Unknown SV type for SV at " << chr << ":" << start << "-" << end << std::endl;
            } else {
                total_count += 1;
            }

            // For complex SVs, split the SV into multiple types (SVTYPE +
            // SVTYPE2)
            SVType sv_type2 = SVType::UNKNOWN;
            if (sv_type == SVType::INV_DEL) {
                sv_type = SVType::DEL;
                sv_type2 = SVType::INV;
            } else if (sv_type == SVType::INV_DUP) {
                sv_type = SVType::DUP;
                sv_type2 = SVType::INV;
            }

            // Deletion
            if (sv_type == SVType::DEL) {
                // Get the deleted sequence from the reference genome, also including the preceding base
                int64_t preceding_pos = (int64_t) std::max(1, (int) start-1);  // Make sure the position is not negative
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
                // Update the position to the preceding base
                int64_t preceding_pos = (int64_t) std::max(1, (int) start-1);  // Make sure the position is not negative
                ref_allele = ref_genome.query(chr, preceding_pos, preceding_pos);
                start = preceding_pos;

                // Update the end position to the same base for duplications and insertions
                if (sv_type == SVType::DUP || sv_type == SVType::INS) {
                    end = start;
                }

                if (sv_type == SVType::INS) {
                    if (alt_allele != "<INS>") {
                        // Insert the reference allele before the insertion
                        alt_allele.insert(0, ref_allele);
                    }
                }
            }

            // Fix ambiguous bases in the reference allele
            const std::string amb_bases = "RYKMSWBDHV";  // Ambiguous bases
            std::bitset<256> amb_bases_bitset;
            for (char base : amb_bases) {
                amb_bases_bitset.set(base);
                amb_bases_bitset.set(std::tolower(base));
            }
            for (char& base : ref_allele) {
                if (amb_bases_bitset.test(base)) {
                    base = 'N';
                }
            }

            // Create the VCF parameter strings
            std::string sv_type_str = getSVTypeString(sv_type);
            std::string sv_type2_str = ".";
            if (sv_type2 != SVType::UNKNOWN) {
                sv_type2_str = getSVTypeString(sv_type2);
            }
            std::string info_str = "END=" + std::to_string(end) + ";SVTYPE=" + sv_type_str + \
                ";SVTYPE2=" + sv_type2_str + ";SVLEN=" + std::to_string(sv_length) + ";SVMETHOD=" + sv_method + ";ALN=" + data_type_str + \
                ";HMM=" + std::to_string(hmm_likelihood) + ";SUPPORT=" + std::to_string(support) + ";CLUSTER=" + std::to_string(cluster_size);

            // std::string info_str = "END=" + std::to_string(end) + ";SVTYPE=" + sv_type_str + 
            //     ";SVLEN=" + std::to_string(sv_length) + ";SVMETHOD=" + sv_method + ";ALN=" + data_type_str + 
            //     ";HMM=" + std::to_string(hmm_likelihood) + ";SUPPORT=" + std::to_string(support) + ";CLUSTER=" + std::to_string(cluster_size);
                
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
    output_vcf += ".gz";
    std::cout << "VCF file created: " << output_vcf << std::endl;
    std::cout << "Index file created: " << output_vcf + ".tbi" << std::endl;

    // Print the number of SV calls skipped
    std::cout << "Finished writing VCF file. Total records: " << total_count << std::endl;
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

    // UPDATE: Only use the start position for the read depth calculation
    // try {
    //     // printMessage("Read depth at end: " + std::to_string(pos_depth_map.at(end)) + " for SV at " + std::to_string(start) + "-" + std::to_string(end) + " with length " + std::to_string(end-start));
    //     read_depth += pos_depth_map.at(end);
    // } catch (const std::out_of_range& e) {
    //     printError("Error: End position " + std::to_string(end) + " not found in depth map.");
    //     // std::cerr << "Warning: End position " << end << " not found in depth map of size " << pos_depth_map.size() << "." << std::endl;
    // }
    // printMessage("Read depth for SV at " + std::to_string(start) + "-" + std::to_string(end) + " with length " + std::to_string(end-start) + ": " + std::to_string(read_depth));
    return read_depth;
}

bool SVCaller::regionOverlaps(const GenomicRegion &a, const GenomicRegion &b)
{
    return a.tid == b.tid && a.start <= b.end && b.start <= a.end;
}
