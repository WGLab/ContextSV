//
// swig_interface.h:
// Declare the C++ functions that will be wrapped by SWIG
//

#ifndef SWIG_INTERFACE_H
#define SWIG_INTERFACE_H

#include "common.h"

#include <string>

int run(std::string bam, std::string snps, std::string outdir, std::string region);

#endif // SWIG_INTERFACE_H
