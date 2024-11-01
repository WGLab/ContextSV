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
{
    this->input_data = &input_data;
}

// Entry point
int ContextSV::run()
{
    // Start the program's timer
    auto start_sv = std::chrono::high_resolution_clock::now();

    // Get the reference genome
    FASTAQuery ref_genome = this->input_data->getRefGenome();

    // Call SVs from long read alignments:
    std::cout << "Running alignment-based SV calling..." << std::endl;
    SVCaller sv_caller(*this->input_data);
    SVData sv_calls = sv_caller.run();

    // Print the total number of SVs called
    // std::cout << "Total SVs called: " << sv_calls.totalCalls() << std::endl;

    // Write SV calls to file
    std::string output_dir = this->input_data->getOutputDir();
    std::cout << "Writing SV calls to file " << output_dir << "/output.vcf..." << std::endl;
    sv_calls.saveToVCF(ref_genome, output_dir);

    // Format and print the time taken to call SVs
    // auto end_sv = std::chrono::high_resolution_clock::now();
    // std::string elapsed_time = getElapsedTime(start_sv, end_sv);
    std::cout << "SV calling complete." << std::endl;
    // int sv_count = sv_calls.totalCalls();
    // std::cout << "Found " << sv_count << " total SVs." << std::endl;
    // std::cout << "SV calling complete. Found " << sv_calls.totalCalls() << " total SVs. Time taken (h:m:s) = " << elapsed_time << std::endl;

    return 0;
}
