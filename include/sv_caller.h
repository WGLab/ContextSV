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
        Common common;

    public:
        SVCaller(Common common);

        // Detect SVs and predict SV type from long read alignments and CNV calls
        SVMap run(CNVMap cnv_calls);
};

#endif // SV_CALLER_H
