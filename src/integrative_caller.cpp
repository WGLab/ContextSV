
#include "integrative_caller.h"
#include "cnv_caller.h"
#include "common.h"

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
    CNVCaller cnv_obj(this->common);
    std::map<int, int> state_sequence =
    cnv_obj.run();

    return 0;
}
