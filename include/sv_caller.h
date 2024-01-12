// SVCaller: Detect SVs from long read alignments

#ifndef SV_CALLER_H
#define SV_CALLER_H

#include "input_data.h"
#include "cnv_data.h"
#include "sv_data.h"

#include <htslib/sam.h>

/// @cond
#include <mutex>
/// @endcond

// SV candidate alignment data (chr, start, end, sequence)
using AlignmentData   = std::tuple<std::string, int64_t, int64_t, std::string>;
using AlignmentVector = std::vector<AlignmentData>;

// Query map (query name, alignment vector)
using QueryMap = std::map<std::string, AlignmentVector>;

class SVCaller {
    private:
        int max_indel_dist = 10;  // Maximum distance between two indels to be considered as a single SV
        int min_sv_size = 50;       // Minimum SV size to be considered
        int min_mapq = 20;          // Minimum mapping quality to be considered
        InputData* input_data;

        // Detect SVs from long read alignments in the CIGAR string
        void detectSVsFromCIGAR(bam_hdr_t* header, bam1_t* alignment, SVData& sv_calls, std::mutex& mtx);

        // Detect SVs from split-read alignments (primary and supplementary)
        SVData detectSVsFromSplitReads(SVData& sv_calls);

        // Detect SVs at a region from long read alignments. This is used for
        // whole genome analysis running in parallel.
        void detectSVsFromRegion(std::string region, SVData &sv_calls, samFile *fp_in, bam_hdr_t *bamHdr, hts_idx_t *idx, std::mutex &callset_mtx, std::mutex &bam_mtx, std::mutex &print_mtx);

        // Read the next alignment from the BAM file in a thread-safe manner
        int readNextAlignment(samFile *fp_in, hts_itr_t *itr, bam1_t *bam1, std::mutex &mtx_bam);

    public:
        SVCaller(InputData& input_data);

        // Detect SVs and predict SV type from long read alignments and CNV calls
        void run(SVData& sv_calls);
};

#endif // SV_CALLER_H
