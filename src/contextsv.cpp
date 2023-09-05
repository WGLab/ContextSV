
#include "contextsv.h"
#include "cnv_caller.h"
#include "sv_caller.h"

#include <htslib/sam.h>
#include <iostream>
#include <string>
#include <vector>


ContextSV::ContextSV(InputData& input_data)
{
    this->input_data = &input_data;
}


// Entry point
int ContextSV::run()
{
    // Call CNVs using the SNP positions
    std::cout << "Calling CNVs..." << std::endl;
    // Get the input data object
    CNVCaller cnv_obj(*this->input_data);
    std::cout << "Starting CNVCaller" << std::endl;
    CNVData cnv_calls = cnv_obj.run();

    // Call SVs from long read alignments and CNV calls
    // Return a map of SV type by start and end position
    // candidate = [chromosome, SV start position], Value = [SV end position, SV type]
    std::cout << "Calling SVs..." << std::endl;
    SVCaller sv_obj(*this->input_data);
    SVData sv_calls = sv_obj.run();

    // Check if the ref genome was set
    if (sv_calls.getRefGenome() == "") {
        std::cout << "Warning: Reference genome not set in SVData" << std::endl;
    }

    // Classify SVs based on CNV calls
    std::cout << "Labeling CNVs..." << std::endl;
    this->labelCNVs(cnv_calls, sv_calls);

    // Write SV calls to file
    std::cout << "Writing SV calls to file..." << std::endl;
    // Get the reference genome
    FASTAQuery ref_genome = this->input_data->getRefGenome();
    std::string output_dir = this->input_data->getOutputDir();
    sv_calls.saveToVCF(ref_genome, output_dir);

    std::cout << "Done!" << std::endl;

    // Integrate CNV and SV calls
    // std::cout << "Integrating CNV and SV calls..." << std::endl;
    // SVMap integrated_calls = integrateCNVs(cnv_calls, sv_calls);

    //filterCNVs(state_sequence, sv_calls);
    //filterCNVs(state_sequence);

    return 0;
}

// Label SVs based on CNV calls
void ContextSV::labelCNVs(CNVData cnv_calls, SVData& sv_calls)
{
    // Iterate over SV calls
    for (auto const& sv_call : sv_calls) {

        SVCandidate candidate = sv_call.first;
        
        // Get the SV coordinates
        std::string chr = std::get<0>(candidate);
        int start_pos = std::get<1>(candidate);
        int end_pos = std::get<2>(candidate);

        // Get the SV type
        SVInfo sv_info = sv_call.second;
        int sv_type = sv_info.first;

        // Get CNV calls within the SV coordinate range and identify the most
        // common call
        int cnv_call = cnv_calls.getMostCommonCNV(chr, start_pos, end_pos);

        // Update the SV call's type
        sv_calls.updateSVType(candidate, cnv_call);

        // Print the updated SV call if the type was changed, and if it is not unknown
        if (cnv_call != -1 && cnv_call != sv_type) {
            std::cout << "Updated SV call from " << sv_type << " to " << cnv_call << std::endl;
        }
    }
}