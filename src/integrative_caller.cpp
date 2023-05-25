
#include "integrative_caller.h"
#include "snv_caller.h"
#include "cnv_caller.h"
#include "common.h"

#include <htslib/sam.h>
#include <iostream>
#include <string>


IntegrativeCaller::IntegrativeCaller(Common common)
{
    this->common = common;
}

/// Entry point
int IntegrativeCaller::run()
{
    // Call SNVs
    //SNVCaller snv_obj;
    //snv_obj.run(filepath);

    // Call CNVs
    CNVCaller cnv_obj(this->common);
    cnv_obj.run();

    return 0;
}
