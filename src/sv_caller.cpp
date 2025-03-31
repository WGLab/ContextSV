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

int SVCaller::readNextAlignment(samFile *fp_in, hts_itr_t *itr, bam1_t *bam1)
{
    std::shared_lock<std::shared_mutex> lock(this->shared_mutex);
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
    }
    bam_hdr_destroy(bamHdr);
    sam_close(fp_in);
    return chromosomes;
}

void SVCaller::findSplitSVSignatures(std::unordered_map<std::string, std::vector<SVCall>> &sv_calls, const InputData &input_data, const std::unordered_map<std::string, std::vector<uint32_t>>& chr_pos_depth_map, const ReferenceGenome& ref_genome)
{
    // Open the BAM file
    std::string bam_filepath = input_data.getLongReadBam();
    samFile *fp_in = sam_open(bam_filepath.c_str(), "r");
    if (!fp_in) {
        printError("ERROR: failed to open " + bam_filepath);
        return;
    }

    // Set maximum thread count
    int thread_count = input_data.getThreadCount();
    hts_set_threads(fp_in, thread_count);
    printMessage("Using " + std::to_string(thread_count) + " threads for split read analysis");

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

    // Alignment data structures
    std::unordered_map<int, std::unordered_map<std::string, PrimaryAlignment>> primary_map;  // TID-> qname -> primary alignment
    std::unordered_map<std::string, std::vector<SuppAlignment>> supp_map;  // qname -> supplementary alignment

    bam1_t *bam1 = bam_init1();
    if (!bam1) {
        printError("ERROR: failed to initialize BAM record");
        return;
    }
    
    // Set the region to the whole genome, or a user-specified chromosome
    hts_itr_t *itr = nullptr;
    if (input_data.isSingleChr()) {
        std::string chr = input_data.getChromosome();
        itr = sam_itr_querys(idx, bamHdr, chr.c_str());
        if (!itr) {
            bam_destroy1(bam1);
            printError("ERROR: failed to create iterator for " + chr);
            return;
        }
    } else {
        itr = sam_itr_queryi(idx, HTS_IDX_START, 0, 0);
        if (!itr) {
            bam_destroy1(bam1);
            printError("ERROR: failed to create iterator for the whole genome");
            return;
        }
    }

    uint32_t primary_count = 0;
    uint32_t supplementary_count = 0;

    // Main loop to process the alignments
    printMessage("Processing alignments from " + bam_filepath);
    uint32_t num_alignments = 0;
    std::unordered_set<int> alignment_tids;  // All unique chromosome IDs
    std::unordered_set<std::string> supp_qnames;  // All unique query names
    while (readNextAlignment(fp_in, itr, bam1) >= 0) {

        // Skip secondary and unmapped alignments, duplicates, QC failures, and low mapping quality
        if (bam1->core.flag & BAM_FSECONDARY || bam1->core.flag & BAM_FUNMAP || bam1->core.flag & BAM_FDUP || bam1->core.flag & BAM_FQCFAIL || bam1->core.qual < this->min_mapq) {
            continue;
        }
        const std::string qname = bam_get_qname(bam1);  // Query template name

        // Process primary alignments
        if (!(bam1->core.flag & BAM_FSUPPLEMENTARY)) {
            // Store chromosome (TID), start, and end positions (1-based) of the
            // primary alignment, and the strand (true for forward, false for
            // reverse)
            std::pair<int, int> qpos = getAlignmentReadPositions(bam1);

            primary_map[bam1->core.tid][qname] = PrimaryAlignment{bam1->core.pos + 1, bam_endpos(bam1), qpos.first, qpos.second, !(bam1->core.flag & BAM_FREVERSE), 0};
            alignment_tids.insert(bam1->core.tid);
            primary_count++;

        // Process supplementary alignments
        } else if (bam1->core.flag & BAM_FSUPPLEMENTARY) {
            // Get the mismatch rate for the read
            const std::string supp_chr = bamHdr->target_name[bam1->core.tid];
            double mismatch_rate = getReadMismatchRate(bam1, supp_chr, ref_genome);

            // Store chromosome (TID), start, and end positions (1-based) of the
            // supplementary alignment, and the strand (true for forward, false
            // for reverse)
            std::pair<int, int> qpos = getAlignmentReadPositions(bam1);
            supp_map[qname].push_back(SuppAlignment{bam1->core.tid, bam1->core.pos + 1, bam_endpos(bam1), qpos.first, qpos.second, !(bam1->core.flag & BAM_FREVERSE), mismatch_rate});
            alignment_tids.insert(bam1->core.tid);
            supp_qnames.insert(qname);
            supplementary_count++;
        }
        num_alignments++;

        if (num_alignments % 1000000 == 0) {
            printMessage("Processed " + std::to_string(num_alignments) + " alignments");
        }
    }

    // Clean up the iterator and alignment
    hts_itr_destroy(itr);
    bam_destroy1(bam1);
    
    // Remove primary alignments without supplementary alignments
    std::unordered_map<int, std::unordered_set<std::string>> to_remove;
    for (auto& chr_primary : primary_map) {
        std::unordered_set<std::string> qnames;
        for (const auto& entry : chr_primary.second) {
            if (supp_qnames.find(entry.first) == supp_qnames.end()) {
                to_remove[chr_primary.first].insert(entry.first);
            }
        }
    }

    int total_removed = 0;
    for (auto& chr_primary : primary_map) {
        // Remove the qnames from the primary map
        total_removed += to_remove[chr_primary.first].size();
        for (const auto& qname : to_remove[chr_primary.first]) {
            chr_primary.second.erase(qname);
        }
    }
    printMessage("Removed " + std::to_string(total_removed) + " primary alignments without supplementary alignments");

    // Process the primary alignments and find SVs
    for (const auto& chr_primary : primary_map) {
        int primary_tid = chr_primary.first;
        std::string chr_name = bamHdr->target_name[primary_tid];
        printMessage("Processing chromosome " + chr_name + " with " + std::to_string(chr_primary.second.size()) + " primary alignments");

        std::vector<SVCall> chr_sv_calls;
        const std::unordered_map<std::string, PrimaryAlignment>& chr_primary_map = chr_primary.second;

        // Identify overlapping primary alignments and cluster endpoints
        std::unique_ptr<IntervalNode> root = nullptr;
        for (const auto& entry : chr_primary_map) {
            const std::string& qname = entry.first;
            const PrimaryAlignment& region = entry.second;
            insert(root, region, qname);
        }

        std::vector<std::vector<std::string>> primary_clusters;
        std::set<std::string> processed;
        for (const auto& entry : chr_primary_map) {
            const std::string& qname = entry.first;
            if (processed.find(qname) != processed.end()) {
                continue;  // Skip already processed primary alignments
            }
            const PrimaryAlignment& primary_aln = entry.second;
            std::vector<std::string> overlap_group;
            findOverlaps(root, primary_aln, overlap_group);
            for (const std::string& qname : overlap_group) {
                processed.insert(qname);
            }
            if (overlap_group.size() > 1) {
                primary_clusters.push_back(overlap_group);
            }
        }

        // For each primary alignment cluster the supplementary alignment start and
        // end positions, keeping the median of the largest cluster
        int current_group = 0;
        int min_length = 2000;
        int max_length = 1000000;
        for (const auto& primary_cluster : primary_clusters) {
            // Determine if the primary alignments are mostly on opposite strands to
            // the corresponding supplementary alignments (potential inversions)
            bool inversion = false;
            int num_primary = (int) primary_cluster.size();
            int num_supp_opposite_strand = 0;
            for (const std::string& qname : primary_cluster) {
                const std::vector<SuppAlignment>& supp_alns = supp_map[qname];
                bool primary_strand = chr_primary_map.at(qname).strand;
                bool has_opposite_strand = false;
                for (const SuppAlignment& supp_aln : supp_alns) {
                    // Analyze if on the same chromosome
                    if (supp_aln.tid == primary_tid && supp_aln.strand != primary_strand) {
                        has_opposite_strand = true;
                    }
                }
                if (has_opposite_strand) {
                    num_supp_opposite_strand++;
                }
            }
            if (static_cast<double>(num_supp_opposite_strand) / static_cast<double>(num_primary) > 0.5) {
                inversion = true;
            }

            // Use DBSCAN to cluster primary alignment start, end positions
            DBSCAN1D dbscan(100, 5);
            current_group++;
            std::vector<int> starts;
            std::vector<int> ends;
            std::vector<bool> primary_strands;
            for (const std::string& qname : primary_cluster) {
                const PrimaryAlignment& primary_aln = chr_primary_map.at(qname);
                starts.push_back(primary_aln.start);
                ends.push_back(primary_aln.end);
                primary_strands.push_back(primary_aln.strand);
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

            // Get the supplementary alignment positions, and also the distances
            // between the primary and supplementary alignments on the read
            std::vector<int> supp_starts;
            std::vector<int> supp_ends;
            std::vector<bool> supp_strands;
            std::vector<int> split_distances;
            std::vector<double> supp_mismatch_rates;
            for (const std::string& qname : primary_cluster) {
                const PrimaryAlignment& primary_aln = chr_primary_map.at(qname);
                const std::vector<SuppAlignment>& supp_alns = supp_map.at(qname);
                for (const SuppAlignment& supp_aln : supp_alns) {
                    if (supp_aln.tid == primary_tid) {
                        // Same chromosome
                        int distance = 0;
                        supp_starts.push_back(supp_aln.start);
                        supp_ends.push_back(supp_aln.end);
                        supp_strands.push_back(supp_aln.strand);
                        supp_mismatch_rates.push_back(supp_aln.mismatch_rate);

                        // Calculate the distance between the primary and supplementary
                        // alignments on the read if on the same chromosome and same
                        // strand
                        if (supp_aln.strand == primary_aln.strand) {
                            // Same strand
                            // Calculate distance (negative if overlapping)
                            if (primary_aln.query_start <= supp_aln.query_start) {
                                distance = supp_aln.query_start - primary_aln.query_end;
                            } else {
                                distance = primary_aln.query_start - supp_aln.query_end;
                            }
                            split_distances.push_back(distance);
                        }

                    } else {
                        // TODO: TRANSLOCATIONS
                    }
                }
            }
            double mean_supp_mismatch_rate = 0.0;
            for (double rate : supp_mismatch_rates) {
                mean_supp_mismatch_rate += rate;
            }
            mean_supp_mismatch_rate /= (double)supp_mismatch_rates.size();

            // Get the largest cluster of supplementary alignment start positions
            dbscan.fit(supp_starts);
            std::vector<int> supp_start_cluster = dbscan.getLargestCluster(supp_starts);

            // Get the largest cluster of supplementary alignment end positions
            dbscan.fit(supp_ends);
            std::vector<int> supp_end_cluster = dbscan.getLargestCluster(supp_ends);

            // Get the largest cluster of split distances
            dbscan.fit(split_distances);
            std::vector<int> split_distance_cluster = dbscan.getLargestCluster(split_distances);

            // Continue if no clusters were found
            if (supp_start_cluster.empty() && supp_end_cluster.empty() && split_distance_cluster.empty()) {
                continue;
            }

            // Use the median of the largest cluster of primary and supplementary
            // alignment start, end positions as the final genome coordinates of the
            // SV
            int primary_pos = -1;
            int primary_pos2 = -1;
            int primary_cluster_size = 0;
            if (primary_start_cluster.size() > primary_end_cluster.size()) {
                std::sort(primary_start_cluster.begin(), primary_start_cluster.end());
                primary_pos = primary_start_cluster[primary_start_cluster.size() / 2];
                primary_cluster_size = primary_start_cluster.size();
            } else if (primary_end_cluster.size() > primary_start_cluster.size()) {
                std::sort(primary_end_cluster.begin(), primary_end_cluster.end());
                primary_pos = primary_end_cluster[primary_end_cluster.size() / 2];
                primary_cluster_size = primary_end_cluster.size();
            } else {
                // Use both positions
                std::sort(primary_start_cluster.begin(), primary_start_cluster.end());
                std::sort(primary_end_cluster.begin(), primary_end_cluster.end());
                primary_pos = primary_start_cluster[primary_start_cluster.size() / 2];
                primary_pos2 = primary_end_cluster[primary_end_cluster.size() / 2];
                primary_cluster_size = primary_start_cluster.size();
            }

            // -------------------------------
            // SPLIT INSERTION DETECTION
            int read_distance = 0;
            if (!split_distance_cluster.empty()) {
                // Use the median of the largest cluster of split distances as the
                // insertion size
                std::sort(split_distance_cluster.begin(), split_distance_cluster.end());
                read_distance = split_distance_cluster[split_distance_cluster.size() / 2];

                // Add an insertion SV call at the primary position
                if (primary_pos != -1 && read_distance > 2000) {
                    if (primary_pos2 != -1) {
                        // If two positions were found, use the 5'most position
                        primary_pos = std::min(primary_pos, primary_pos2);
                    }
                    //int read_depth = this->getReadDepth(chr_pos_depth_map.at(chr_name), primary_pos);
                    SVCall sv_candidate(primary_pos, primary_pos + (read_distance-1), SVType::INS, getSVTypeSymbol(SVType::INS), SVDataType::SPLITDIST1, Genotype::UNKNOWN, 0.0, 0, mean_supp_mismatch_rate, primary_cluster_size);
                    
                    // Print if end position = 162908547
                    if (primary_pos + (read_distance - 1) == 162908547) {
                        printMessage("[TEST] Adding insertion SV candidate at " + chr_name + ":" + std::to_string(primary_pos) + "-" + std::to_string(primary_pos + (read_distance - 1)) + " with length " + std::to_string(read_distance));
                    }
                    addSVCall(chr_sv_calls, sv_candidate);
                }
            }

            // --------------------------------

            // Get the supplementary alignment positions
            int supp_pos = -1;
            int supp_pos2 = -1;
            int supp_cluster_size = 0;
            int supp_best_start = -1;
            int supp_best_end = -1;
            if (!supp_start_cluster.empty()) {
                std::sort(supp_start_cluster.begin(), supp_start_cluster.end());
                supp_best_start = supp_start_cluster[supp_start_cluster.size() / 2];
            }
            if (!supp_end_cluster.empty()) {
                std::sort(supp_end_cluster.begin(), supp_end_cluster.end());
                supp_best_end = supp_end_cluster[supp_end_cluster.size() / 2];
            }

            if (supp_start_cluster.size() > supp_end_cluster.size()) {
                supp_pos = supp_best_start;
                supp_cluster_size = supp_start_cluster.size();
            } else if (supp_end_cluster.size() > supp_start_cluster.size()) {
                supp_pos = supp_best_end;
                supp_cluster_size = supp_end_cluster.size();
            } else if (supp_best_end == -1 && supp_best_start == -1) {
                // Use both positions. This has been shown to occur in some nested SVs
                supp_pos = supp_best_start;
                supp_pos2 = supp_best_end;
                supp_cluster_size = supp_start_cluster.size();
            }

            // Store the inversion as the supplementary start and end positions
            if (supp_best_start != -1 && supp_best_end != -1) {
                if (inversion && std::abs(supp_best_start - supp_best_end) >= 50) {
                    //int read_depth = this->getReadDepth(chr_pos_depth_map.at(chr_name), std::min(supp_best_start, supp_best_end));

                    // Print if end position = 162908547
                    if (std::max(supp_best_start, supp_best_end) == 162908547) {
                        printMessage("[TEST] Adding inversion SV candidate at " + chr_name + ":" + std::to_string(std::min(supp_best_start, supp_best_end)) + "-" + std::to_string(std::max(supp_best_start, supp_best_end)) + " with length " + std::to_string(std::abs(supp_best_start - supp_best_end)));
                    }

                    SVCall sv_candidate(std::min(supp_best_start, supp_best_end), std::max(supp_best_start, supp_best_end), SVType::INV, getSVTypeSymbol(SVType::INV), SVDataType::SUPPINV, Genotype::UNKNOWN, 0.0, 0, mean_supp_mismatch_rate, supp_cluster_size);
                    addSVCall(chr_sv_calls, sv_candidate);
                }
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
            int cluster_size = std::max(primary_cluster_size, supp_cluster_size);

            // If the read distance is < 30bp while the SV is > 2kb, then this is a
            // potential deletion
            if (std::abs(read_distance) < 30 && sv_length > 2000 && sv_length <= 1000000) {

                // Add an inversion call if necessary
                //int read_depth = this->getReadDepth(chr_pos_depth_map.at(chr_name), sv_start);
                if (inversion) {
                    SVCall sv_candidate(sv_start, sv_end, SVType::INV, getSVTypeSymbol(SVType::INV), SVDataType::SPLITINV, Genotype::UNKNOWN, 0.0, 0, mean_supp_mismatch_rate, cluster_size);
                    addSVCall(chr_sv_calls, sv_candidate);
                } else {
                    SVCall sv_candidate(sv_start, sv_end, SVType::DEL, getSVTypeSymbol(SVType::DEL), SVDataType::SPLITDIST2, Genotype::UNKNOWN, 0.0, 0, mean_supp_mismatch_rate, cluster_size);
                    addSVCall(chr_sv_calls, sv_candidate);
                }
            }

            // Add a dummy SV call for CNV detection
            else if (sv_length >= min_length && sv_length <= max_length) {
                SVType sv_type = inversion ? SVType::INV : SVType::UNKNOWN;
                std::string alt = (sv_type == SVType::INV) ? "<INV>" : ".";
                //int read_depth = this->getReadDepth(chr_pos_depth_map.at(chr_name), sv_start);
                SVCall sv_candidate(sv_start, sv_end, sv_type, alt, SVDataType::SPLIT, Genotype::UNKNOWN, 0.0, 0, mean_supp_mismatch_rate, cluster_size);
                addSVCall(chr_sv_calls, sv_candidate);
            }
        }

        // Combine SVs with identical start and end positions, and sum the cluster
        // sizes
        printMessage("Combining SVs with identical start positions");
        std::sort(chr_sv_calls.begin(), chr_sv_calls.end(), [](const SVCall& a, const SVCall& b) {
            return a.start < b.start || (a.start == b.start && a.end < b.end);
        });
        
        // Merge duplicate SV calls with identical start positions
        mergeDuplicateSVs(chr_sv_calls);

        printMessage("Merged SVs:");
        for (const auto& sv : chr_sv_calls) {
            printMessage(" - " + getSVTypeSymbol(sv.sv_type) + " at " + chr_name + ":" + std::to_string(sv.start) + "-" + std::to_string(sv.end) + " with length " + std::to_string(sv.end - sv.start + 1) + " and cluster size " + std::to_string(sv.cluster_size));
        }

        sv_calls[chr_name] = std::move(chr_sv_calls);

        // Print the number of merged SV calls
        printMessage(chr_name + ": Found " + std::to_string(sv_calls[chr_name].size()) + " SV candidates");
    }
}

void SVCaller::findCIGARSVs(samFile* fp_in, hts_idx_t* idx, bam_hdr_t* bamHdr, const std::string& region, std::vector<SVCall>& sv_calls, const std::vector<uint32_t>& pos_depth_map)
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

        // Skip secondary and unmapped alignments, duplicates, QC failures, and
        // low mapping quality, and supplementary alignments
        if (bam1->core.flag & BAM_FSECONDARY || bam1->core.flag & BAM_FUNMAP || bam1->core.flag & BAM_FDUP || bam1->core.flag & BAM_FQCFAIL || bam1->core.qual < this->min_mapq || bam1->core.flag & BAM_FSUPPLEMENTARY) {
            continue;
        }

        // Process the alignment
        // bool primary = !(bam1->core.flag & BAM_FSUPPLEMENTARY);
        this->processCIGARRecord(bamHdr, bam1, sv_calls, pos_depth_map);
    }

    // Clean up the iterator and alignment
    hts_itr_destroy(itr);
    bam_destroy1(bam1);
}

double SVCaller::getReadMismatchRate(bam1_t *alignment, const std::string& chr, const ReferenceGenome & ref_genome)
{
    uint32_t* cigar = bam_get_cigar(alignment);  // CIGAR array
    int cigar_len = alignment->core.n_cigar;
    uint32_t query_pos = 0;
    uint32_t pos = (uint32_t)alignment->core.pos;
    uint32_t aln_start = pos;
    uint32_t end = (uint32_t)bam_endpos(alignment) - 1;  // Rightmost position of the alignment in the reference genome (0-based)

    // Get the reference sequence
    std::string_view ref_seq = ref_genome.query(chr, pos + 1, end + 1);

    // Loop through the CIGAR string and calculate the number of matches and
    // mismatches
    int match_count = 0;
    int mismatch_count = 0;
    for (int i = 0; i < cigar_len; i++) {
        int op_len = bam_cigar_oplen(cigar[i]);  // CIGAR operation length
        int op = bam_cigar_op(cigar[i]);  // CIGAR operation
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            for (int j = 0; j < op_len; j++) {
                char base = seq_nt16_str[bam_seqi(bam_get_seq(alignment), query_pos + j)];
                if (base == ref_seq[pos - aln_start + j]) {
                    match_count++;
                } else {
                    mismatch_count++;
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

    // Calculate the mismatch rate
    double mismatch_rate = 0.0;
    if (match_count + mismatch_count > 0) {
        mismatch_rate = static_cast<double>(mismatch_count) / static_cast<double>(match_count + mismatch_count);
    }
    return mismatch_rate;
}

void SVCaller::processCIGARRecord(bam_hdr_t *header, bam1_t *alignment, std::vector<SVCall> &sv_calls, const std::vector<uint32_t> &pos_depth_map)
{
    std::string chr = header->target_name[alignment->core.tid];  // Chromosome name
    uint32_t aln_start = (uint32_t)alignment->core.pos;  // Leftmost position of the alignment in the reference genome (0-based)
    uint32_t pos = aln_start;
    // uint32_t end = (uint32_t)bam_endpos(alignment) - 1;  // Rightmost position of the alignment in the reference genome (0-based)

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

    std::vector<SVCall> cigar_sv_calls;
    cigar_sv_calls.reserve(1000);
    for (int i = 0; i < cigar_len; i++) {
        int op_len = bam_cigar_oplen(cigar[i]);  // CIGAR operation length
        int op = bam_cigar_op(cigar[i]);  // CIGAR operation
        if (op_len >= 50) {
            
            // Process the CIGAR operation
            if (op == BAM_CINS) {

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

                // Add as an insertion
                uint32_t ins_pos = pos + 1;
                uint32_t ins_end = ins_pos + op_len - 1;
                //int read_depth = this->getReadDepth(pos_depth_map, ins_pos);
                
                // Determine the ALT allele format based on small vs. large insertion
                std::string alt_allele = "<INS>";
                if (op_len <= 50) {
                    alt_allele = ins_seq_str;
                }
                SVCall sv_call(ins_pos, ins_end, SVType::INS, alt_allele, SVDataType::CIGARINS, Genotype::UNKNOWN, default_lh, 0, 1, 0);
                cigar_sv_calls.emplace_back(sv_call);
            
            // Process clipped bases as potential insertions
            } else if (op == BAM_CSOFT_CLIP) {
                // Soft-clipped bases are considered as potential insertions
                // Skip if the position exceeds the reference genome length
                if (pos + 1 >= pos_depth_map.size()) {
                    continue;
                }

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
                
                // Add as an insertion
                uint32_t ins_pos = pos + 1;
                uint32_t ins_end = ins_pos + op_len - 1;
                //int read_depth = this->getReadDepth(pos_depth_map, ins_pos);

                // Determine the ALT allele format based on small vs. large insertion
                std::string alt_allele = "<INS>";
                if (op_len <= 50) {
                    alt_allele = ins_seq_str;
                }
                SVCall sv_call(ins_pos, ins_end, SVType::INS, alt_allele, SVDataType::CIGARCLIP, Genotype::UNKNOWN, default_lh, 0, 0.0, 0);
                cigar_sv_calls.emplace_back(sv_call);  // Commented for testing

            // Check if the CIGAR operation is a deletion
            } else if (op == BAM_CDEL) {

                ref_pos = pos+1;
                ref_end = ref_pos + op_len -1;
                //int read_depth = this->getReadDepth(pos_depth_map, ref_pos);
                // SVCall sv_call(ref_pos, ref_end, SVType::DEL, "<DEL>",
                // "CIGARDEL", "./.", default_lh, read_depth, 1, 0);
                SVCall sv_call(ref_pos, ref_end, SVType::DEL, getSVTypeSymbol(SVType::DEL), SVDataType::CIGARDEL, Genotype::UNKNOWN, default_lh, 0, 1, 0);
                // addSVCall(sv_calls, sv_call);
                // addSVCall(cigar_sv_calls, sv_call);
                cigar_sv_calls.emplace_back(sv_call);
            }
            
            // For matches, calculate the sequence identity
            // } else if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            //     if (ref_seq.size() < static_cast<size_t>(op_len)) {
            //         printError("ERROR: reference sequence length is less than the CIGAR operation length");
            //         continue;
            //     }

            //     // printMessage("Calculating sequence identity for matches");
            //     for (int j = 0; j < op_len; j++) {
            //         char base = seq_nt16_str[bam_seqi(bam_get_seq(alignment), query_pos + j)];
            //         if (base == ref_seq[pos - aln_start + j]) {
            //             match_count++;
            //         } else {
            //             mismatch_count++;
            //         }
            //     }
            // }
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

    // Get the read name
    // std::string read_name = bam_get_qname(alignment);

    // If read name starts with c08844a5 then print the read name and the number
    // of matches and mismatches
    // if (read_name.find("c08844a5") != std::string::npos) {
    //     printMessage(read_name + ": matches=" + std::to_string(match_count) + ", mismatches=" + std::to_string(mismatch_count) + ", mismatches/length=" + std::to_string((double)mismatch_count / (double)(match_count + mismatch_count)));
    // }
    // double mismatch_rate = (double)mismatch_count / (double)(match_count + mismatch_count);
    // if (mismatch_rate > 0) {
    //     printMessage("Read name: " + read_name + ", mismatch rate: " + std::to_string(mismatch_rate) + ", matches: " + std::to_string(match_count) + ", mismatches: " + std::to_string(mismatch_count));
    // }
    // read_mismatch_rates[read_name] = mismatch_rate;
    // printMessage("Completed processing read: " + read_name);

    // Set the mismatch rate for all SVs from this read, and add the SV calls
    for (SVCall& sv_call : cigar_sv_calls) {
        // sv_call.mismatch_rate = mismatch_rate;
        addSVCall(sv_calls, sv_call);
    }
}

std::pair<int, int> SVCaller::getAlignmentReadPositions(bam1_t *alignment)
{
    int query_start = -1;
    int query_end = 0;
    uint32_t* cigar = bam_get_cigar(alignment);
    int cigar_len = alignment->core.n_cigar;
    for (int i = 0; i < cigar_len; i++) {
        int op_len = bam_cigar_oplen(cigar[i]);
        int op = bam_cigar_op(cigar[i]);

        // Set the query start position to the first non-soft clip operation
        if (query_start == -1 && (op == BAM_CMATCH || op == BAM_CINS || op == BAM_CEQUAL || op == BAM_CDIFF)) {
            query_start = query_end;  // First valid query position
        }
        
        // https://github.com/samtools/htslib/blob/develop/htslib/sam.h:
        // bam_cigar_type(o) (BAM_CIGAR_TYPE>>((o)<<1)&3) // bit 1: consume query; bit 2: consume reference
        if (op == BAM_CMATCH || op == BAM_CINS || op == BAM_CSOFT_CLIP || op == BAM_CEQUAL || op == BAM_CDIFF) {
            query_end += op_len;
        }
    }

    if (query_start == -1) {
        query_start = 0;
    }

    return std::make_pair(query_start, query_end);
}

void SVCaller::processChromosome(const std::string& chr, std::vector<SVCall>& chr_sv_calls, const InputData& input_data, const std::vector<uint32_t>& chr_pos_depth_map, double mean_chr_cov)
{
    // Open the BAM file
    std::string bam_filepath = input_data.getLongReadBam();
    samFile *fp_in = sam_open(bam_filepath.c_str(), "r");
    if (!fp_in) {
        printError("ERROR: failed to open " + bam_filepath);
        return;
    }
    hts_set_threads(fp_in, 1);

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

    // Get DBSCAN parameters
    double dbscan_epsilon = input_data.getDBSCAN_Epsilon();
    int dbscan_min_pts = 5;
    double dbscan_min_pts_pct = input_data.getDBSCAN_MinPtsPct();
    if (dbscan_min_pts_pct > 0.0) {
        dbscan_min_pts = (int)std::ceil(mean_chr_cov * dbscan_min_pts_pct);
        printMessage(chr + ": Mean chr. cov.: " + std::to_string(mean_chr_cov) + " (DBSCAN min. pts.= " + std::to_string(dbscan_min_pts) + ", min. pts. pct.= " + std::to_string(dbscan_min_pts_pct) + ")");
    } 

    // -----------------------------------------------------------------------
    // Detect SVs from the CIGAR strings
    printMessage(chr + ": CIGAR SVs...");
    this->findCIGARSVs(fp_in, idx, bamHdr, chr, chr_sv_calls, chr_pos_depth_map);

    printMessage(chr + ": Merging CIGAR...");
    mergeSVs(chr_sv_calls, dbscan_epsilon, dbscan_min_pts, false);

    int region_sv_count = getSVCount(chr_sv_calls);
    printMessage(chr + ": Found " + std::to_string(region_sv_count) + " SV candidates in the CIGAR string");
}

void SVCaller::run(const InputData& input_data)
{
    bool cigar_svs = true;
    bool split_svs = true;

    // Set up the reference genome
    printMessage("Loading the reference genome...");
    const std::string ref_filepath = input_data.getRefGenome();
    std::shared_mutex ref_mutex;  // Dummy mutex (remove later)
    ReferenceGenome ref_genome(ref_mutex);
    ref_genome.setFilepath(ref_filepath);

    // Get the chromosomes
    std::vector<std::string> chromosomes;
    if (input_data.isSingleChr()) {
        // Get the chromosome from the user input argument
        chromosomes.push_back(input_data.getChromosome());
    } else {
        // Get the chromosomes from the input BAM file
        chromosomes = this->getChromosomes(input_data.getLongReadBam());
    }
    
    // Read the HMM from the file
    std::string hmm_filepath = input_data.getHMMFilepath();
    std::cout << "Reading HMM from file: " << hmm_filepath << std::endl;
    const CHMM& hmm = ReadCHMM(hmm_filepath.c_str());

    // Set up the JSON output file for CNV data
    const std::string& json_fp = input_data.getCNVOutputFile();

    // Calculate the mean chromosome coverage and generate the position depth
    // maps for each chromosome (I/O is multi-threaded, which is more efficient
    // than per-chromosome multi-threading in this case)
    std::shared_mutex shared_mutex;
    CNVCaller cnv_caller(shared_mutex);
    std::unordered_map<std::string, std::vector<uint32_t>> chr_pos_depth_map;
    std::unordered_map<std::string, double> chr_mean_cov_map;
    const std::string bam_filepath = input_data.getLongReadBam();
    int chr_thread_count = input_data.getThreadCount();

    // Initialize the chromosome position depth map and mean coverage map
    for (const auto& chr : chromosomes) {
        uint32_t chr_len = ref_genome.getChromosomeLength(chr);
        if (chr_len == 0) {
            printError("Chromosome " + chr + " not found in reference genome");
            continue;
        }
        chr_pos_depth_map[chr] = std::vector<uint32_t>(chr_len+1, 0);  // 1-based index
        chr_mean_cov_map[chr] = 0.0;
    }
    cnv_caller.calculateMeanChromosomeCoverage(chromosomes, chr_pos_depth_map, chr_mean_cov_map, bam_filepath, chr_thread_count);

    // Remove chromosomes with no reads (mean coverage is zero)
    std::vector<std::string> null_chr;
    for (const auto& chr : chromosomes) {
        if (chr_mean_cov_map[chr] == 0.0) {
            null_chr.push_back(chr);
        }
    }
    for (const auto& chr : null_chr) {
        printMessage("Removing chromosome " + chr + " with no reads...");
        chromosomes.erase(std::remove(chromosomes.begin(), chromosomes.end(), chr), chromosomes.end());
    }
    std::unordered_map<std::string, std::vector<SVCall>> whole_genome_sv_calls;
    int current_chr = 0;
    int total_chr_count = chromosomes.size();

    if (cigar_svs) {
        // Use multi-threading across chromosomes. If a single chromosome is
        // specified, use a single main thread (multi-threading is used for file I/O)
        int thread_count = 1;
        if (!input_data.isSingleChr()) {
            thread_count = input_data.getThreadCount();
            std::cout << "Using " << thread_count << " threads for chr processing..." << std::endl;
        }
        ThreadPool pool(thread_count);
        auto process_chr = [&](const std::string& chr) {
            try {
                std::vector<SVCall> sv_calls;
                InputData chr_input_data = input_data;  // Use a thread-local copy
                this->processChromosome(chr, sv_calls, chr_input_data, chr_pos_depth_map[chr], chr_mean_cov_map[chr]);
                {
                    std::shared_lock<std::shared_mutex> lock(this->shared_mutex);
                    whole_genome_sv_calls[chr] = std::move(sv_calls);
                }
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
                process_chr(chr);
            }));
        }

        // Wait for all tasks to complete
        for (auto& future : futures) {
            try {
                current_chr++;
                future.get();
            } catch (const std::exception& e) {
                printError("Error processing chromosome task: " + std::string(e.what()));
            } catch (...) {
                printError("Unknown error processing chromosome task.");
            }
        }
        printMessage("All tasks have finished.");

        // -------------------------------------------------------
        // Run copy number variant predictions on the SVs detected from the
        // CIGAR string, using a minimum CNV length threshold
        current_chr = 0;
        printMessage("Running copy number predictions on CIGAR SVs...");
        for (auto& entry : whole_genome_sv_calls) {
            current_chr++;
            const std::string& chr = entry.first;
            std::vector<SVCall>& sv_calls = entry.second;
            if (sv_calls.size() > 0) {
                printMessage("(" + std::to_string(current_chr) + "/" + std::to_string(total_chr_count) + ") Running copy number predictions on " + chr + "...");
                cnv_caller.runCIGARCopyNumberPrediction(chr, sv_calls, hmm, chr_mean_cov_map[chr], chr_pos_depth_map[chr], input_data);
            }
        }
        // -------------------------------------------------------
    }
    
    if (split_svs) {
        // Identify split-SV signatures
        printMessage("Identifying split-SV signatures...");
        std::unordered_map<std::string, std::vector<SVCall>> whole_genome_split_sv_calls;
        this->findSplitSVSignatures(whole_genome_split_sv_calls, input_data, chr_pos_depth_map, ref_genome);

        printMessage("Running copy number predictions on split-read SVs...");
        current_chr = 0;
        for (auto& entry : whole_genome_split_sv_calls) {
            const std::string& chr = entry.first;
            std::vector<SVCall>& sv_calls = entry.second;

            if (sv_calls.size() > 0) {
                current_chr++;
                printMessage("(" + std::to_string(current_chr) + "/" + std::to_string(total_chr_count) + ") Running copy number predictions on " + chr + " with " + std::to_string(sv_calls.size()) + " SV candidates...");
                this->runSplitReadCopyNumberPredictions(chr, sv_calls, cnv_caller, hmm, chr_mean_cov_map[chr], chr_pos_depth_map[chr], input_data);
            }
        }

        printMessage("Merging split-read SVs...");
        int min_pts = 2;
        for (auto& entry : whole_genome_split_sv_calls) {
            std::vector<SVCall>& sv_calls = entry.second;
            mergeSVs(sv_calls, input_data.getDBSCAN_Epsilon(), min_pts, true);
        }

        printMessage("Unifying SVs...");
        for (auto& entry : whole_genome_split_sv_calls) {
            const std::string& chr = entry.first;
            std::vector<SVCall>& sv_calls = entry.second;
            whole_genome_sv_calls[chr].insert(whole_genome_sv_calls[chr].end(), sv_calls.begin(), sv_calls.end());
        }
    }

    if (input_data.getSaveCNVData()) {
        closeJSON(json_fp);
    }

    // Print the total number of SVs detected for each chromosome
    uint32_t total_sv_count = 0;
    for (const auto& entry : whole_genome_sv_calls) {
        std::string chr = entry.first;
        int sv_count = getSVCount(entry.second);
        total_sv_count += sv_count;
        printMessage("Total SVs detected for " + chr + ": " + std::to_string(sv_count));
    }
    printMessage("Total SVs detected: " + std::to_string(total_sv_count));

    // Save to VCF
    std::cout << "Saving SVs to VCF..." << std::endl;
    const std::string output_dir = input_data.getOutputDir();
    this->saveToVCF(whole_genome_sv_calls, output_dir, ref_genome, chr_pos_depth_map);
}

void SVCaller::findOverlaps(const std::unique_ptr<IntervalNode> &root, const PrimaryAlignment &query, std::vector<std::string> &result)
{
    if (!root) return;

    // If overlapping, add to result
    if (query.start <= root->region.end && query.end >= root->region.start)
        result.push_back(root->qname);

    // If left subtree may have overlaps, search left
    if (root->left && root->left->max_end >= query.start)
        findOverlaps(root->left, query, result);

    // Always check the right subtree
    findOverlaps(root->right, query, result);
}

void SVCaller::insert(std::unique_ptr<IntervalNode> &root, const PrimaryAlignment &region, std::string qname)
{
    if (!root) {
        root = std::make_unique<IntervalNode>(region, qname);
        return;
    }

    if (region.start < root->region.start)
    {
        insert(root->left, region, qname);
    } else {
        insert(root->right, region, qname);
    }

    // Update max_end
    root->max_end = std::max(root->max_end, region.end);
}

// Run copy number predictions on the SVs detected from the split reads
void SVCaller::runSplitReadCopyNumberPredictions(const std::string& chr, std::vector<SVCall>& split_sv_calls, const CNVCaller& cnv_caller, const CHMM& hmm, double mean_chr_cov, const std::vector<uint32_t>& pos_depth_map, const InputData& input_data)
{
    std::vector<SVCall> additional_calls;
    for (auto& sv_candidate : split_sv_calls) {
        std::tuple<double, SVType, Genotype, bool> result = cnv_caller.runCopyNumberPrediction(chr, hmm, sv_candidate.start, sv_candidate.end, mean_chr_cov, pos_depth_map, input_data);
        double supp_lh = std::get<0>(result);
        SVType supp_type = std::get<1>(result);
        Genotype genotype = std::get<2>(result);

        // For inversions with copy-neutral support, update the HMM likelihood
        if (supp_type == SVType::NEUTRAL && sv_candidate.sv_type == SVType::INV) {
            sv_candidate.hmm_likelihood = supp_lh;
        }

        // Update the SV type if the support is not neutral or unknown
        else if (supp_type != SVType::UNKNOWN && supp_type != SVType::NEUTRAL) {
            // Update information if the SV call is unknown
            if (sv_candidate.sv_type == SVType::UNKNOWN) {
                sv_candidate.sv_type = supp_type;
                sv_candidate.alt_allele = getSVTypeSymbol(supp_type);  // Update the ALT allele format
                sv_candidate.data_type = SVDataType::HMM;
                sv_candidate.genotype = genotype;
                sv_candidate.hmm_likelihood = supp_lh;

                // Print if end position = 162908547
                if (sv_candidate.end == 162908547) {
                    printMessage("SV at " + chr + ":" + std::to_string(sv_candidate.start) + "-" + std::to_string(sv_candidate.end) +
                        " updated to type " + getSVTypeSymbol(supp_type) +
                        " with likelihood " + std::to_string(supp_lh) +
                        " and genotype " + getGenotypeString(genotype));
                }
            // Add an additional SV call if the type is different
            } else if (sv_candidate.sv_type != supp_type) {
                SVCall new_sv_call = sv_candidate;  // Copy the original SV call
                new_sv_call.sv_type = supp_type;
                new_sv_call.alt_allele = getSVTypeSymbol(supp_type);  // Update the ALT allele format
                new_sv_call.data_type = SVDataType::HMM;
                new_sv_call.genotype = genotype;
                new_sv_call.hmm_likelihood = supp_lh;

                // Print if end position = 162908547
                if (new_sv_call.end == 162908547) {
                    printMessage("Additional SV at " + chr + ":" + std::to_string(new_sv_call.start) + "-" + std::to_string(new_sv_call.end) +
                        " with type " + getSVTypeSymbol(supp_type) +
                        " and likelihood " + std::to_string(supp_lh) +
                        " and genotype " + getGenotypeString(genotype));
                }
                additional_calls.push_back(new_sv_call);
            }
        }
    }

    // Add the additional SV calls to the original list, replacing any existing
    // ones
    for (auto& new_sv_call : additional_calls) {
        bool found = false;
        for (auto& existing_sv_call : split_sv_calls) {
            if (existing_sv_call.start == new_sv_call.start && existing_sv_call.end == new_sv_call.end &&
                existing_sv_call.sv_type == new_sv_call.sv_type) {
                // Update the existing SV call with the new one
                existing_sv_call = new_sv_call;
                found = true;
                break;
            }
        }
        if (!found) {
            addSVCall(split_sv_calls, new_sv_call);  // Add as a new SV call
        }
    }
}

void SVCaller::saveToVCF(const std::unordered_map<std::string, std::vector<SVCall>>& sv_calls, const std::string& output_dir, const ReferenceGenome& ref_genome, const std::unordered_map<std::string, std::vector<uint32_t>>& chr_pos_depth_map) const
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
    std::cout << "Formatting VCF header..." << std::endl;
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
        "##INFO=<ID=MISMATCH,Number=1,Type=Float,Description=\"Mismatch rate\">",
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
    int unclassified_svs = 0;
    int filtered_svs = 0;
    for (const auto& pair : sv_calls) {
        std::string chr = pair.first;
        const std::vector<SVCall>& sv_calls = pair.second;
        std::cout << "Saving SV calls for " << chr << "..." << std::endl;
        for (const auto& sv_call : sv_calls) {
            // Get the SV candidate and SV info
            uint32_t start = sv_call.start;
            uint32_t end = sv_call.end;
            SVType sv_type = sv_call.sv_type;
            // std::string genotype = sv_call.genotype;
            // std::string data_type_str = sv_call.data_type;
            // std::string alt_allele = sv_call.alt_allele;
            std::string genotype = getGenotypeString(sv_call.genotype);
            std::string data_type_str = getSVDataTypeString(sv_call.data_type);
            std::string alt_allele = sv_call.alt_allele;
            double hmm_likelihood = sv_call.hmm_likelihood;
            int sv_length = end - start + 1;
            int cluster_size = sv_call.cluster_size;
            //int read_depth = sv_call.read_depth;
            std::string ref_allele = ".";
            double mismatch_rate = sv_call.mismatch_rate;
            std::string filter = "PASS";

            // If the SV type is unknown, print a warning and skip
            if (sv_type == SVType::UNKNOWN || sv_type == SVType::NEUTRAL) {
                unclassified_svs += 1; 
                continue;
            } else {
                total_count += 1;
            }

            // Deletion
            if (sv_type == SVType::DEL) {
                // Get the deleted sequence from the reference genome, also including the preceding base
                uint32_t preceding_pos = (uint32_t) std::max(1, (int) start-1);  // Make sure the position is not negative
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

                if (sv_type == SVType::INS) {
                    // Update the position to the preceding base
                    int64_t preceding_pos = (int64_t) std::max(1, (int) start-1);  // Make sure the position is not negative
                    ref_allele = ref_genome.query(chr, preceding_pos, preceding_pos);
                    start = preceding_pos;

                    if (alt_allele != "<INS>") {
                        // Insert the reference allele before the insertion
                        alt_allele.insert(0, ref_allele);
                    }
                    end = start;  // Update the end position to the same base

                } else {
                    ref_allele = "N";  // Convention for INV and DUP
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
            
            // Get read depth
            int read_depth = this->getReadDepth(chr_pos_depth_map.at(chr), start);

            // Create the VCF parameter strings
            std::string sv_type_str = getSVTypeString(sv_type);
            // std::string info_str = "END=" + std::to_string(end) + ";SVTYPE=" + sv_type_str + ";SVLEN=" + std::to_string(sv_length) + ";SVMETHOD=" + sv_method + ";ALN=" + data_type_str + ";HMM=" + std::to_string(hmm_likelihood) + ";SUPPORT=" + std::to_string(read_depth) + ";CLUSTER=" + std::to_string(cluster_size);                
            std::string info_str = "END=" + std::to_string(end) + ";SVTYPE=" + sv_type_str + ";SVLEN=" + std::to_string(sv_length) + ";SVMETHOD=" + sv_method + ";ALN=" + data_type_str + ";HMM=" + std::to_string(hmm_likelihood) + ";SUPPORT=" + std::to_string(read_depth) + ";CLUSTER=" + std::to_string(cluster_size) + ";MISMATCH=" + std::to_string(mismatch_rate);
            std::string format_str = "GT:DP";
            std::string sample_str = genotype + ":" + std::to_string(read_depth);
            std::vector<std::string> samples = {sample_str};

            // Write the SV call to the file (CHROM, POS, ID, REF, ALT, QUAL,
            // FILTER, INFO, FORMAT, SAMPLES)
            vcf_stream << chr << "\t" << start << "\t" << "." << "\t" << ref_allele << "\t" << alt_allele << "\t" << "." << "\t" << filter << "\t" << info_str << "\t" << format_str << "\t" << samples[0] << std::endl;
            // vcf_stream << chr << "\t" << start << "\t" << "." << "\t" << ref_allele << "\t" << alt_allele << "\t" << "." << "\t" << "PASS" << "\t" << info_str << "\t" << format_str << "\t" << samples[0] << std::endl;
        }
    }
    vcf_stream.close();
    std::cout << "Saved SV calls to " << output_vcf << std::endl;

    // Print the number of SV calls skipped
    std::cout << "Finished writing VCF file. Total records: " << total_count << std::endl;
    if (unclassified_svs > 0) {
        std::cout << "Total unclassified SVs: " << unclassified_svs << std::endl;
    }
    printMessage("Total PASS filtered SVs: " + std::to_string(filtered_svs));
}

int SVCaller::getReadDepth(const std::vector<uint32_t>& pos_depth_map, uint32_t start) const
{
    int read_depth = 0;
    try {
        read_depth += pos_depth_map.at(start);
    } catch (const std::out_of_range& e) {
        // Occurs with clipped reads (insertion evidence) that are outside the
        // range of the depth map
        printError("Warning: Read depth for position " + std::to_string(start) + " is out of range of size " + std::to_string(pos_depth_map.size()));
    }

    return read_depth;
}
