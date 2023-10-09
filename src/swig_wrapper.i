/*
SWIG wrapper for C++ code.
*/

%module contextsv

// Include header
%{
#include "swig_interface.h"
#include "input_data.h"
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

// Expose the InputData class
%include "input_data.h"

// Include functions
int run(InputData input_data);

