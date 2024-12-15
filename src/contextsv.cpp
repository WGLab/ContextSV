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

ContextSV::ContextSV(InputData& input_data)
    : input_data(input_data)  // Initialize the input data
{
}

int ContextSV::run()
{
    SVCaller sv_caller(this->input_data); 
    sv_caller.run();

    return 0;
}
