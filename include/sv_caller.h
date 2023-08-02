//
// sv_caller.cpp:
// Detect SVs from long read alignments
//

#ifndef SV_CALLER_H
#define SV_CALLER_H

#include "common.h"

class SVCaller {
    private:
        Common common;

    public:
        SVCaller(Common common);

        // Detect SVs and return SV type by start and end position
        std::map<std::pair<char *, int>, std::pair<int, int>> run();
};

#endif // SV_CALLER_H
