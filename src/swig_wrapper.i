/*
SWIG wrapper for C++ code.
*/

%module contextsv

// Include header
%{
#include "cli.h"
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
int run(const std::string& bam, const std::string& snps, const std::string& outdir, const std::string& region);
