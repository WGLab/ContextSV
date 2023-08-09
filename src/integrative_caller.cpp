
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
    CNVMap cnv_calls = cnv_obj.run();

    // Call SVs from long read alignments and CNV calls
    // Return a map of SV type by start and end position
    // Key = [chromosome, SV start position], Value = [SV end position, SV type]
    std::cout << "Calling SVs..." << std::endl;
    SVCaller sv_obj(this->common);
    SVMap sv_calls = sv_obj.run(cnv_calls);

    // Integrate CNV and SV calls
    // std::cout << "Integrating CNV and SV calls..." << std::endl;
    // SVMap integrated_calls = integrateCNVs(cnv_calls, sv_calls);

    //filterCNVs(state_sequence, sv_calls);
    //filterCNVs(state_sequence);

    return 0;
}
