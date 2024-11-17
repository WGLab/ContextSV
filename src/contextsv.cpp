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
    FASTAQuery ref_genome = this->input_data->getRefGenome();  // Load the reference genome
    SVCaller sv_caller(*this->input_data);  // Create an SV caller object
    SVData sv_calls = sv_caller.run();  // Run the SV caller
    std::string output_dir = this->input_data->getOutputDir();  // Get the output directory
    
    std::cout << "Writing SV calls to file " << output_dir << "/output.vcf..." << std::endl;
    sv_calls.saveToVCF(ref_genome, output_dir);  // Save the SV calls to a VCF file
    std::cout << "SV calling complete." << std::endl;

    return 0;
}
