//
// swig_interface.h:
// Declare the C++ functions that will be wrapped by SWIG
//

#ifndef SWIG_INTERFACE_H
#define SWIG_INTERFACE_H

#include "input_data.h"

/// @cond
#include <string>
/// @endcond

int run(std::string bam_fp, std::string ref_fp, std::string snps_fp, std::string outdir, std::string region, std::string chr_cov);

#endif // SWIG_INTERFACE_H
