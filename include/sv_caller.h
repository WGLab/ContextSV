// SVCaller: Detect SVs from long read alignments

#ifndef SV_CALLER_H
#define SV_CALLER_H

#include "cnv_caller.h"
#include "input_data.h"
#include "sv_object.h"
#include "fasta_query.h"

#include <htslib/sam.h>

/// @cond
// #include <mutex>
#include <shared_mutex>
#include <unordered_map>
#include <future>
/// @endcond

struct GenomicRegion {
    int tid;
    hts_pos_t start;
    hts_pos_t end;
    int query_start;
    int query_end;
    bool strand;
    int cluster_size;  // Number of alignments used for this region
};

struct PrimaryAlignment {
    hts_pos_t start;
    hts_pos_t end;
    int query_start;
    int query_end;
    bool strand;
    int cluster_size;  // Number of alignments used for this region
};

struct SuppAlignment {
    int tid;
    hts_pos_t start;
    hts_pos_t end;
    int query_start;
    int query_end;
    bool strand;
    int cluster_size;  // Number of alignments used for this region
};

struct SplitSignature {
    int tid;
    hts_pos_t start;
    hts_pos_t end;
    bool strand;
    hts_pos_t query_start;
    hts_pos_t query_end;
};

// Interval Tree Node
struct IntervalNode {
    PrimaryAlignment region;
    std::string qname;
    hts_pos_t max_end;  // To optimize queries
    std::unique_ptr<IntervalNode> left;
    std::unique_ptr<IntervalNode> right;

    IntervalNode(PrimaryAlignment r, std::string name)
        : region(r), qname(name), max_end(r.end), left(nullptr), right(nullptr) {}
};

void insert(std::unique_ptr<IntervalNode>& root, const PrimaryAlignment& region, std::string qname) {
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

void findOverlaps(const std::unique_ptr<IntervalNode>& root, const PrimaryAlignment& query, std::vector<std::string>& result) {
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

struct MismatchData {
    uint32_t query_start;
    uint32_t query_end;
    std::vector<int> match_map;
};

class SVCaller {
    private:
        int min_mapq = 20;          // Minimum mapping quality to be considered
        mutable std::shared_mutex shared_mutex;  // Shared mutex for thread safety

        std::vector<std::string> getChromosomes(const std::string& bam_filepath);

        // void findSplitCNVBreakpoints(samFile* fp_in, hts_idx_t* idx, bam_hdr_t* bamHdr, const std::string& region, std::vector<SVCall>& sv_calls);

        void findSplitSVSignatures(std::unordered_map<std::string, std::vector<SVCall>>& sv_calls, const InputData& input_data);

        void findSplitReadSVs(std::unordered_map<std::string, std::vector<SVCall>>& sv_calls, const ReferenceGenome& ref_genome, const InputData& input_data);

        // Process a single CIGAR record and find candidate SVs
        void processCIGARRecord(bam_hdr_t* header, bam1_t* alignment, std::vector<SVCall>& sv_calls, bool is_primary, const std::vector<uint32_t>& pos_depth_map, const ReferenceGenome& ref_genome);

        std::pair<int, int> getAlignmentReadPositions(bam1_t* alignment);

        void processChromosome(const std::string& chr, const CHMM& hmm, std::vector<SVCall>& combined_sv_calls, const InputData& input_data, const ReferenceGenome& ref_genome, const std::vector<uint32_t>& chr_pos_depth_map, double mean_chr_cov, std::vector<SVCall>& split_sv_calls);

        // Detect SVs at a region from long read alignments. This is used for
        // whole genome analysis running in parallel.
        // RegionData detectSVsFromRegion(std::string region);
        void findCIGARSVs(samFile* fp_in, hts_idx_t* idx, bam_hdr_t* bamHdr, const std::string& region, std::vector<SVCall>& sv_calls, const std::vector<uint32_t>& pos_depth_map, const ReferenceGenome& ref_genome);
 
        // Read the next alignment from the BAM file in a thread-safe manner
        int readNextAlignment(samFile *fp_in, hts_itr_t *itr, bam1_t *bam1);

        // Detect SVs from split alignments
        // void detectSVsFromSplitReads(const std::string& region, samFile* fp_in, hts_idx_t* idx, bam_hdr_t* bamHdr, std::vector<SVCall>& split_sv_calls, const CNVCaller& cnv_caller, const CHMM& hmm, double mean_chr_cov, const std::vector<uint32_t>& pos_depth_map, const InputData& input_data);

        // Calculate the mismatch rate given a map of query positions to
        // match/mismatch (1/0) values within a specified range of the query
        // sequence
        double calculateMismatchRate(const MismatchData& mismatch_data);

        void runSplitReadCopyNumberPredictions(const std::string& chr, std::vector<SVCall>& split_sv_calls, const CNVCaller &cnv_caller, const CHMM &hmm, double mean_chr_cov, const std::vector<uint32_t> &pos_depth_map, const InputData &input_data);

        void saveToVCF(const std::unordered_map<std::string, std::vector<SVCall>> &sv_calls, const std::string &output_dir, const ReferenceGenome &ref_genome) const;

        // Calculate the read depth (INFO/DP) for a region
        int calculateReadDepth(const std::vector<uint32_t>& pos_depth_map, uint32_t start, uint32_t end);

    public:
        // Constructor with no arguments
        SVCaller() = default;

        // Detect SVs and predict SV type from long read alignments and CNV calls
        void run(const InputData& input_data);
};

#endif // SV_CALLER_H
