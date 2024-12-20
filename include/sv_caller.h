// SVCaller: Detect SVs from long read alignments

#ifndef SV_CALLER_H
#define SV_CALLER_H

#include "cnv_caller.h"
#include "input_data.h"
#include "sv_object.h"
#include "fasta_query.h"

#include <htslib/sam.h>

/// @cond
#include <mutex>
#include <unordered_map>
#include <future>
/// @endcond

struct GenomicRegion {
    int tid;
    hts_pos_t start;
    hts_pos_t end;
    bool strand;
};

struct MismatchData {
    uint32_t query_start;
    uint32_t query_end;
    std::vector<int> match_map;
};

class SVCaller {
    private:
        int min_sv_size = 50;       // Minimum SV size to be considered
        int min_mapq = 20;          // Minimum mapping quality to be considered

        // void getAlignmentMismatchMap(samFile* fp_in, hts_idx_t* idx, bam_hdr_t* bamHdr, const GenomicRegion& region, MismatchData& mismatch_data, bool is_primary) const;

        void getSplitAlignments(samFile* fp_in, hts_idx_t* idx, bam_hdr_t* bamHdr, const std::string& region, std::unordered_map<std::string, GenomicRegion>& primary_map, std::unordered_map<std::string, std::vector<GenomicRegion>>& supp_map) const;

        // Detect SVs from the CIGAR string of a read alignment, and return the
        // mismatch rate, and the start and end positions of the query sequence
        void detectSVsFromCIGAR(bam_hdr_t* header, bam1_t* alignment, std::vector<SVCall>& sv_calls, bool is_primary, const std::vector<uint32_t>& pos_depth_map, const ReferenceGenome& ref_genome) const;

        void processChromosome(const std::string& chr, const CHMM& hmm, std::vector<SVCall>& combined_sv_calls, const InputData& input_data, const ReferenceGenome& ref_genome, std::mutex& snp_mutex, std::mutex& pfb_mutex);

        // Detect SVs at a region from long read alignments. This is used for
        // whole genome analysis running in parallel.
        // RegionData detectSVsFromRegion(std::string region);
        void detectCIGARSVs(samFile* fp_in, hts_idx_t* idx, bam_hdr_t* bamHdr, const std::string& region, std::vector<SVCall>& sv_calls, const std::vector<uint32_t>& pos_depth_map, const ReferenceGenome& ref_genome) const;
 
        // Read the next alignment from the BAM file in a thread-safe manner
        int readNextAlignment(samFile *fp_in, hts_itr_t *itr, bam1_t *bam1) const;

        // Detect SVs from split alignments
        void detectSVsFromSplitReads(const std::string& region, samFile* fp_in, hts_idx_t* idx, bam_hdr_t* bamHdr, std::vector<SVCall>& sv_calls, const CNVCaller& cnv_caller, const CHMM& hmm, double mean_chr_cov, const std::vector<uint32_t>& pos_depth_map, const InputData& input_data, std::mutex& snp_mutex, std::mutex& pfb_mutex) const;

        // Calculate the mismatch rate given a map of query positions to
        // match/mismatch (1/0) values within a specified range of the query
        // sequence
        double calculateMismatchRate(const MismatchData& mismatch_data) const;

        void saveToVCF(const std::unordered_map<std::string, std::vector<SVCall>>& sv_calls, const std::string& output_dir, const ReferenceGenome& ref_genome) const;

        void trimOverlappingAlignments(GenomicRegion& primary_alignment, GenomicRegion& supp_alignment, const MismatchData& primary_mismatches, const MismatchData& supp_mismatches) const;

        // Calculate the read depth (INFO/DP) for a region
        int calculateReadDepth(const std::vector<uint32_t>& pos_depth_map, uint32_t start, uint32_t end) const;

    public:
        // Constructor with no arguments
        SVCaller() = default;

        // Detect SVs and predict SV type from long read alignments and CNV calls
        void run(const InputData& input_data);
};

#endif // SV_CALLER_H
