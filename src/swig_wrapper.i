/*
SWIG wrapper for C++ code.
*/

%module contextsv

// Include header
%{
#include "swig_interface.h"
%}

// Set up types
%include "std_string.i"
%include "stdint.i"

// Define an output handler to redirect stdout to Python
%feature("python:output") {
    // Redirect stdout to Python
    void print_stdout(const std::string& msg) {
        PySys_WriteStdout("%s", msg.c_str());
    }
}

// Include functions
int run(    \
    const std::string& bam_fp,  \
    const std::string& ref_fp,  \
    const std::string& snps_fp, \
    const std::string& outdir,  \
    const std::string& region,  \
    const std::string& chr_cov, \
    const std::string& pfb_fp,  \
    int thread_count    \
    );
