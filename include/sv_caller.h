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
        //int max_indel_dist = 1000;  // Maximum distance between two indels to
        //be considered as a single SV
        int max_indel_dist = 10;  // Maximum distance between two indels to be considered as a single SV
        //int min_sv_size = 50;       // Minimum SV size to be considered
        int min_sv_size = 30;       // Minimum SV size to be considered
        Common common;

    public:
        SVCaller(Common common);

        // Detect SVs and predict SV type from long read alignments and CNV calls
        SVMap run(CNVMap cnv_calls);

        // Detect SVs from long read alignments
        SVMap detectSVs();

        // Merge SVs
        SVMap mergeSVs(SVMap sv_calls);

        // Predict SV type from CNV calls
        SVMap predictSVType(SVMap sv_calls, CNVMap cnv_calls);
};

#endif // SV_CALLER_H
