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

void SVCaller::findSplitSVSignatures(std::unordered_map<std::string, std::vector<SVCall>> &sv_calls, const InputData &input_data)
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
            // Store chromosome (TID), start, and end positions (1-based) of the
            // supplementary alignment, and the strand (true for forward, false
            // for reverse)
            std::pair<int, int> qpos = getAlignmentReadPositions(bam1);
            supp_map[qname].push_back(SuppAlignment{bam1->core.tid, bam1->core.pos + 1, bam_endpos(bam1), qpos.first, qpos.second, !(bam1->core.flag & BAM_FREVERSE), 0});
            alignment_tids.insert(bam1->core.tid);
            supp_qnames.insert(qname);
            supplementary_count++;
        }
        num_alignments++;

        if (num_alignments % 1000000 == 0) {
            printMessage("Processed " + std::to_string(num_alignments) + " alignments");
        }
    }

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
                        } else {
                            // TODO: INVERSIONS                       
                        }
                    } else {
                        // TODO: TRANSLOCATIONS
                    }
                }
            }

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
                    SVCall sv_candidate(primary_pos, primary_pos + (read_distance-1), SVType::INS, "<INS>", "SPLITINS", "./.", 0.0, 0, 0, primary_cluster_size);
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
                    SVCall sv_candidate(std::min(supp_best_start, supp_best_end), std::max(supp_best_start, supp_best_end), SVType::INV, "<INV>", "SUPPINV", "./.", 0.0, 0, 0, supp_cluster_size);
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
                SVCall sv_candidate(sv_start, sv_end, SVType::DEL, ".", "SPLITDEL", "./.", 0.0, 0, 0, cluster_size);
                addSVCall(chr_sv_calls, sv_candidate);

                // Add an inversion call if necessary
                if (inversion) {
                    SVCall sv_candidate(sv_start, sv_end, SVType::INV, "<INV>", "INVDEL", "./.", 0.0, 0, 0, cluster_size);
                    addSVCall(chr_sv_calls, sv_candidate);
                }
            }

            // Add a dummy SV call for CNV detection
            else if (sv_length >= min_length && sv_length <= max_length) {
                SVType sv_type = inversion ? SVType::INV : SVType::UNKNOWN;
                std::string alt = (sv_type == SVType::INV) ? "<INV>" : ".";
                SVCall sv_candidate(sv_start, sv_end, sv_type, alt, "PRIMSUPP", "./.", 0.0, 0, 0, cluster_size);
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
        sv_calls[chr_name] = std::move(chr_sv_calls);

        // Print the number of merged SV calls
        printMessage(chr_name + ": Found " + std::to_string(sv_calls[chr_name].size()) + " SV candidates");

        // Print all SV calls
        for (const SVCall& sv_call : sv_calls[chr_name]) {
            printMessage("SV: " + std::to_string(sv_call.start) + "-" + std::to_string(sv_call.end) + " " + getSVTypeString(sv_call.sv_type) + ", length: " + std::to_string(sv_call.end - sv_call.start + 1) + ", cluster size: " + std::to_string(sv_call.cluster_size) + ", group: " + std::to_string(current_group));
        }
    }
}

void SVCaller::findCIGARSVs(samFile* fp_in, hts_idx_t* idx, bam_hdr_t* bamHdr, const std::string& region, std::vector<SVCall>& sv_calls, const std::vector<uint32_t>& pos_depth_map, const ReferenceGenome& ref_genome)
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
        this->processCIGARRecord(bamHdr, bam1, sv_calls, primary, pos_depth_map, ref_genome);
    }

    // Clean up the iterator and alignment
    hts_itr_destroy(itr);
    bam_destroy1(bam1);
}

void SVCaller::processCIGARRecord(bam_hdr_t *header, bam1_t *alignment, std::vector<SVCall> &sv_calls, bool is_primary, const std::vector<uint32_t> &pos_depth_map, const ReferenceGenome &ref_genome)
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
                        int read_depth = this->getReadDepth(pos_depth_map, bp1);
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
                        int read_depth = this->getReadDepth(pos_depth_map, bp1);
                        SVCall sv_call(bp1, bp2, SVType::DUP, "<DUP>", "RSEQSIM", "./.", default_lh, read_depth, 1, 0);
                        addSVCall(sv_calls, sv_call);
                        continue;
                    }
                }

                // Add as an insertion
                uint32_t ins_pos = pos + 1;
                uint32_t ins_end = ins_pos + op_len - 1;
                int read_depth = this->getReadDepth(pos_depth_map, ins_pos-1);
                
                // Determine the ALT allele format based on small vs. large insertion
                std::string alt_allele = "<INS>";
                if (op_len <= 50) {
                    alt_allele = ins_seq_str;
                }
                SVCall sv_call(ins_pos, ins_end, SVType::INS, alt_allele, "CIGARINS", "./.", default_lh, read_depth, 1, 0);
                addSVCall(sv_calls, sv_call);

            // Check if the CIGAR operation is a deletion
            } else if (op == BAM_CDEL && is_primary) {

                ref_pos = pos+1;
                ref_end = ref_pos + op_len -1;
                int read_depth = this->getReadDepth(pos_depth_map, ref_pos);
                SVCall sv_call(ref_pos, ref_end, SVType::DEL, "<DEL>", "CIGARDEL", "./.", default_lh, read_depth, 1, 0);
                addSVCall(sv_calls, sv_call);
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

void SVCaller::processChromosome(const std::string& chr, std::vector<SVCall>& chr_sv_calls, const InputData& input_data, const ReferenceGenome& ref_genome, const std::vector<uint32_t>& chr_pos_depth_map, double mean_chr_cov)
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
    this->findCIGARSVs(fp_in, idx, bamHdr, chr, chr_sv_calls, chr_pos_depth_map, ref_genome);

    printMessage(chr + ": Merging CIGAR...");
    mergeSVs(chr_sv_calls, dbscan_epsilon, dbscan_min_pts);

    int region_sv_count = getSVCount(chr_sv_calls);
    printMessage(chr + ": Found " + std::to_string(region_sv_count) + " SV candidates in the CIGAR string");
}

void SVCaller::run(const InputData& input_data)
{
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
            std::vector<SVCall> split_sv_calls;
            InputData chr_input_data = input_data;  // Use a thread-local copy
            this->processChromosome(chr, sv_calls, chr_input_data, ref_genome, chr_pos_depth_map[chr], chr_mean_cov_map[chr]);
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
            // printMessage("Processing chromosome " + chr);
            process_chr(chr);
        }));
    }

    // // Wait for all tasks to complete
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

    // Identify split-SV signatures
    printMessage("Identifying split-SV signatures...");
    std::unordered_map<std::string, std::vector<SVCall>> whole_genome_split_sv_calls;
    this->findSplitSVSignatures(whole_genome_split_sv_calls, input_data);

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
    
    printMessage("Unifying SVs...");
    for (auto& entry : whole_genome_split_sv_calls) {
        const std::string& chr = entry.first;
        std::vector<SVCall>& sv_calls = entry.second;
        whole_genome_sv_calls[chr].insert(whole_genome_sv_calls[chr].end(), sv_calls.begin(), sv_calls.end());
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
    this->saveToVCF(whole_genome_sv_calls, output_dir, ref_genome);
}

// Run copy number predictions on the SVs detected from the split reads
void SVCaller::runSplitReadCopyNumberPredictions(const std::string& chr, std::vector<SVCall>& split_sv_calls, const CNVCaller& cnv_caller, const CHMM& hmm, double mean_chr_cov, const std::vector<uint32_t>& pos_depth_map, const InputData& input_data)
{
    std::vector<SVCall> processed_calls;
    for (const auto& sv_candidate : split_sv_calls) {
        std::tuple<double, SVType, std::string, bool> result = cnv_caller.runCopyNumberPrediction(chr, hmm, sv_candidate.start, sv_candidate.end, mean_chr_cov, pos_depth_map, input_data);
        double supp_lh = std::get<0>(result);
        SVType supp_type = std::get<1>(result);
        std::string genotype = std::get<2>(result);
        if (supp_type != SVType::UNKNOWN && supp_type != SVType::NEUTRAL) {
            int read_depth = this->getReadDepth(pos_depth_map, sv_candidate.start);
            std::string alt_allele = "<" + getSVTypeString(supp_type) + ">";
            SVCall sv_call(sv_candidate.start, sv_candidate.end, supp_type, alt_allele, "SPLIT", genotype, supp_lh, read_depth, 1, sv_candidate.cluster_size);
            processed_calls.push_back(sv_call);
        }
    }

    // Insert the copy number predictions back into the split SV calls
    printMessage("Inserting CNV calls...");
    split_sv_calls.insert(split_sv_calls.end(), processed_calls.begin(), processed_calls.end());
    mergeDuplicateSVs(split_sv_calls);
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
    std::cout << "Formatting VCF header..." << std::endl;
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
    int unclassified_svs = 0;
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
                unclassified_svs += 1;
                continue;
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
                } else {
                    ref_allele = "N";  // Convention for INV and DUP
                }
                end = start;  // Update the end position to the same base
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
                
            std::string format_str = "GT:DP";
            std::string sample_str = genotype + ":" + std::to_string(read_depth);
            std::vector<std::string> samples = {sample_str};

            // Write the SV call to the file (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLES)
            vcf_stream << chr << "\t" << start << "\t" << "." << "\t" << ref_allele << "\t" << alt_allele << "\t" << "." << "\t" << "PASS" << "\t" << info_str << "\t" << format_str << "\t" << samples[0] << std::endl;
        }
    }
    vcf_stream.close();
    std::cout << "Saved SV calls to " << output_vcf << std::endl;

    // Print the number of SV calls skipped
    std::cout << "Finished writing VCF file. Total records: " << total_count << std::endl;
    if (unclassified_svs > 0) {
        std::cout << "Total unclassified SVs: " << unclassified_svs << std::endl;
    }
}

int SVCaller::getReadDepth(const std::vector<uint32_t>& pos_depth_map, uint32_t start)
{
    int read_depth = 0;
    try {
        read_depth += pos_depth_map.at(start);
    } catch (const std::out_of_range& e) {
        printError("Error: Start position " + std::to_string(start) + " not found in depth map.");
    }

    return read_depth;
}
