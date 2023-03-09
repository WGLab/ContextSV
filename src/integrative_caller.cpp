
#include "integrative_caller.h"
#include "snv_caller.h"
#include "cnv_caller.h"

#include <iostream>
#include <string>


IntegrativeCaller::IntegrativeCaller()
= default;

/// Entry point
int IntegrativeCaller::run(std::string filepath)
{
    // Call SNVs
    SNVCaller snv_obj;
    snv_obj.run(filepath);

    // Call CNVs
    CNVCaller cnv_obj;
    cnv_obj.run(filepath, snv_obj);

    return 0;
}
