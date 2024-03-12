// SVCaller: Detect SVs from long read alignments

#ifndef SV_CALLER_H
#define SV_CALLER_H

#include "input_data.h"
#include "cnv_data.h"
#include "sv_data.h"

#include <htslib/sam.h>

/// @cond
#include <mutex>
#include <unordered_map>
#include <future>
/// @endcond

// SV candidate alignment data (chr, start, end, sequence)
using AlignmentData   = std::tuple<std::string, int64_t, int64_t, std::string>;
using AlignmentVector = std::vector<AlignmentData>;

// Query map (query name, alignment vector)
using PrimaryMap = std::unordered_map<std::string, AlignmentData>;
using SuppMap = std::unordered_map<std::string, AlignmentVector>;
using RegionData = std::tuple<SVData, PrimaryMap, SuppMap>;

class SVCaller {
    private:
        int max_indel_dist = 10;  // Maximum distance between two indels to be considered as a single SV
        int min_sv_size = 50;       // Minimum SV size to be considered
        int min_mapq = 20;          // Minimum mapping quality to be considered
        InputData* input_data;
        std::mutex bam_mtx;  // Mutex for locking the BAM file
        std::mutex print_mtx;  // Mutex for locking printing to stdout
        std::mutex query_mtx;  // Mutex for locking the query map
        std::mutex sv_mtx;  // Mutex for locking the SV data

        // Detect SVs from long read alignments in the CIGAR string
        void detectSVsFromCIGAR(bam_hdr_t* header, bam1_t* alignment, SVData& sv_calls);

        // Detect SVs at a region from long read alignments. This is used for
        // whole genome analysis running in parallel.
        RegionData detectSVsFromRegion(std::string region);
        // SVData detectSVsFromRegion(std::string region, samFile *fp_in, bam_hdr_t *bamHdr, hts_idx_t *idx);
 
        // Read the next alignment from the BAM file in a thread-safe manner
        int readNextAlignment(samFile *fp_in, hts_itr_t *itr, bam1_t *bam1);

        void detectSVsFromSplitReads(SVData& sv_calls, PrimaryMap& primary_map, SuppMap& supp_map);

    public:
        SVCaller(InputData& input_data);

        // Detect SVs and predict SV type from long read alignments and CNV calls
        SVData run();
};

#endif // SV_CALLER_H
