//
// sv_caller.cpp:
// Detect SVs from long read alignments
//

#ifndef SV_CALLER_H
#define SV_CALLER_H

#include "common.h"
#include "sv_map.h"

class SVCaller {
    private:
        Common common;

    public:
        SVCaller(Common common);

        // Detect SVs and return SV type by start and end position
        SVMap run();
};

#endif // SV_CALLER_H
