#include "contextsv.h"
#include "sv_caller.h"

#include <htslib/sam.h>

/// @cond
#include <iostream>
#include <string>
#include <vector>
#include <chrono>

#include "utils.h"
/// @endcond


int ContextSV::run(const InputData& input_data) const
{
    SVCaller sv_caller;
    sv_caller.run(input_data);

    return 0;
}
