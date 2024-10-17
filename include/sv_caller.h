// SVCaller: Detect SVs from long read alignments

#ifndef SV_CALLER_H
#define SV_CALLER_H

#include "cnv_caller.h"
#include "input_data.h"
#include "cnv_data.h"
#include "sv_data.h"

#include <htslib/sam.h>

/// @cond
#include <mutex>
#include <unordered_map>
#include <future>
/// @endcond

// SV candidate alignment data (chr, start, end, sequence, query start, query
// end, mismatch map)
using AlignmentData   = std::tuple<std::string, int64_t, int64_t, std::string, int32_t, int32_t, std::unordered_map<int, int>>;
using AlignmentVector = std::vector<AlignmentData>;

// Query map (query name, alignment vector)
using PrimaryMap = std::unordered_map<std::string, AlignmentData>;
using SuppMap = std::unordered_map<std::string, AlignmentVector>;
using RegionData = std::tuple<SVData, PrimaryMap, SuppMap>;

class SVCaller {
    private:
        int min_sv_size = 50;       // Minimum SV size to be considered
        int min_mapq = 20;          // Minimum mapping quality to be considered
        InputData* input_data;
        std::mutex sv_mtx;  // Mutex for locking the SV data

        // Detect SVs from the CIGAR string of a read alignment, and return the
        // mismatch rate, and the start and end positions of the query sequence
        std::tuple<std::unordered_map<int, int>, int32_t, int32_t> detectSVsFromCIGAR(bam_hdr_t* header, bam1_t* alignment, SVData& sv_calls, bool is_primary);

        // Detect SVs at a region from long read alignments. This is used for
        // whole genome analysis running in parallel.
        RegionData detectSVsFromRegion(std::string region);
 
        // Read the next alignment from the BAM file in a thread-safe manner
        int readNextAlignment(samFile *fp_in, hts_itr_t *itr, bam1_t *bam1);

        // Detect SVs from split alignments
        void detectSVsFromSplitReads(SVData& sv_calls, PrimaryMap& primary_map, SuppMap& supp_map, CNVCaller& cnv_caller);

        // Calculate the mismatch rate given a map of query positions to
        // match/mismatch (1/0) values within a specified range of the query
        // sequence
        double calculateMismatchRate(std::unordered_map<int, int>& mismatch_map, int32_t start, int32_t end);

    public:
        SVCaller(InputData& input_data);

        // Detect SVs and predict SV type from long read alignments and CNV calls
        SVData run();
};

#endif // SV_CALLER_H
