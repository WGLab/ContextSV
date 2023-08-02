
#include "integrative_caller.h"
#include "common.h"
#include "cnv_caller.h"
#include "sv_caller.h"

#include <htslib/sam.h>
#include <iostream>
#include <string>
#include <vector>


IntegrativeCaller::IntegrativeCaller(Common common)
{
    this->common = common;
}

/// Entry point
int IntegrativeCaller::run()
{
    // Call CNVs using the SNP positions
    std::cout << "Calling CNVs..." << std::endl;
    CNVCaller cnv_obj(this->common);
    std::map<int, int> state_sequence =
    cnv_obj.run();

    // Call SVs from long read alignments
    std::cout << "Calling SVs..." << std::endl;
    SVCaller sv_obj(this->common);

    // Return a map of SV type by start and end position (key=[chromosome, SV
    // start position], value=[SV end position, SV type])
    std::map<std::pair<char *, int>, std::pair<int, int>> sv_calls =
    sv_obj.run();

    // Integrate CNV and SV calls
    // TODO
    std::cout << "Integrating CNV and SV calls..." << std::endl;
    //filterCNVs(state_sequence);

    return 0;
}
