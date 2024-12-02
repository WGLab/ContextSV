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

// SV candidate alignment data (chr, start, end, sequence, query start, query
// end, mismatch map, strand)
using AlignmentData   = std::tuple<std::string, uint32_t, uint32_t, uint32_t, uint32_t, std::vector<int>, bool>;
using AlignmentVector = std::vector<AlignmentData>;

// Query map (query name, alignment vector)
using PrimaryMap = std::unordered_map<std::string, AlignmentData>;
using SuppMap = std::unordered_map<std::string, AlignmentVector>;
// using RegionData = std::tuple<SVData, PrimaryMap, SuppMap>;

class SVCaller {
    private:
        int min_sv_size = 50;       // Minimum SV size to be considered
        int min_mapq = 20;          // Minimum mapping quality to be considered
        InputData& input_data;

        // Detect SVs from the CIGAR string of a read alignment, and return the
        // mismatch rate, and the start and end positions of the query sequence
        void detectSVsFromCIGAR(bam_hdr_t* header, bam1_t* alignment, std::vector<SVCall>& sv_calls, std::tuple<std::vector<int>, uint32_t, uint32_t>& query_info, bool is_primary, const std::vector<uint32_t>& pos_depth_map);

        void processChromosome(const std::string& chr, const std::string& bam_filepath, const CHMM& hmm, std::vector<SVCall>& combined_sv_calls, int min_cnv_length);

        // Detect SVs at a region from long read alignments. This is used for
        // whole genome analysis running in parallel.
        // RegionData detectSVsFromRegion(std::string region);
        void detectCIGARSVs(samFile* fp_in, hts_idx_t* idx, bam_hdr_t* bamHdr, const std::string& region, std::vector<SVCall>& sv_calls, PrimaryMap& primary_alignments, SuppMap& supplementary_alignments, const std::vector<uint32_t>& pos_depth_map);
 
        // Read the next alignment from the BAM file in a thread-safe manner
        int readNextAlignment(samFile *fp_in, hts_itr_t *itr, bam1_t *bam1);

        // Detect SVs from split alignments
        void detectSVsFromSplitReads(std::vector<SVCall>& sv_calls, PrimaryMap& primary_map, SuppMap& supp_map, CNVCaller& cnv_caller, const CHMM& hmm, double mean_chr_cov, const std::vector<uint32_t>& pos_depth_map);

        // Calculate the mismatch rate given a map of query positions to
        // match/mismatch (1/0) values within a specified range of the query
        // sequence
        double calculateMismatchRate(const std::vector<int>& mismatch_map, int32_t start, int32_t end);

        void saveToVCF(const std::unordered_map<std::string, std::vector<SVCall>>& sv_calls);

        void trimOverlappingAlignments(AlignmentData& primary_alignment, AlignmentData& supp_alignment);

        // Calculate the read depth (INFO/DP) for a region
        int calculateReadDepth(const std::vector<uint32_t>& pos_depth_map, uint32_t start, uint32_t end);

    public:
        explicit SVCaller(InputData& input_data);

        // Detect SVs and predict SV type from long read alignments and CNV calls
        void run();
};

#endif // SV_CALLER_H
