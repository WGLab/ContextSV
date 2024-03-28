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

// Set up the namespace
%include "input_data.h"

// Include functions
int run(InputData input_data);

