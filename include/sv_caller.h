//
// sv_caller.cpp:
// Detect SVs from long read alignments
//

#ifndef SV_CALLER_H
#define SV_CALLER_H

#include "common.h"
#include "cnv_map.h"
#include "sv_map.h"

class SVCaller {
    private:
        // Aliases
        typedef std::tuple<std::string, int, int> AlignmentLocation;
        typedef std::vector<AlignmentLocation> AlignmentVector;
        typedef std::map<std::string, AlignmentVector> QueryMap;

        //int max_indel_dist = 1000;  // Maximum distance between two indels to
        //be considered as a single SV
        int max_indel_dist = 10;  // Maximum distance between two indels to be considered as a single SV
        //int min_sv_size = 50;       // Minimum SV size to be considered
        //int min_sv_size = 30;       // Minimum SV size to be considered
        int min_sv_size = 50;       // Minimum SV size to be considered
        int min_mapq = 20;          // Minimum mapping quality to be considered
        Common common;

    public:
        SVCaller(Common common);

        // Detect SVs and predict SV type from long read alignments and CNV calls
        SVMap run(CNVMap cnv_calls);

        // Detect SVs from long read alignments in the CIGAR string
        SVMap detectSVsFromCIGAR();

        // Detect SVs from split-read alignments (primary and supplementary)
        SVMap detectSVsFromSplitReads();
};

#endif // SV_CALLER_H
