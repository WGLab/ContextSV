// SVCaller: Detect SVs from long read alignments

#ifndef SV_CALLER_H
#define SV_CALLER_H

#include "input_data.h"
#include "cnv_data.h"
#include "sv_data.h"

#include <htslib/sam.h>

// SV candidate alignment data (chr, start, end, sequence)
using AlignmentData   = std::tuple<std::string, int64_t, int64_t, std::string>;
using AlignmentVector = std::vector<AlignmentData>;

// Query map (query name, alignment vector)
using QueryMap = std::map<std::string, AlignmentVector>;

class SVCaller {
    private:
        //int max_indel_dist = 1000;  // Maximum distance between two indels to
        //be considered as a single SV
        int max_indel_dist = 10;  // Maximum distance between two indels to be considered as a single SV
        //int min_sv_size = 50;       // Minimum SV size to be considered
        //int min_sv_size = 30;       // Minimum SV size to be considered
        int min_sv_size = 50;       // Minimum SV size to be considered
        int min_mapq = 20;          // Minimum mapping quality to be considered
        InputData* input_data;

        // Detect SVs from long read alignments in the CIGAR string
        void detectSVsFromCIGAR(bam_hdr_t* header, bam1_t* alignment, SVData& sv_calls);
        //void detectSVsFromCIGAR(SVData& sv_calls, std::string chr, int32_t pos, uint32_t* cigar, int cigar_len, bool debug_mode);

        // Detect SVs from split-read alignments (primary and supplementary)
        SVData detectSVsFromSplitReads(SVData& sv_calls);

    public:
        SVCaller(InputData& input_data);

        // Detect SVs and predict SV type from long read alignments and CNV calls
        void run(SVData& sv_calls);
};

#endif // SV_CALLER_H
