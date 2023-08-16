
#include "contextsv.h"
#include "common.h"
#include "cnv_caller.h"
#include "sv_caller.h"

#include <htslib/sam.h>
#include <iostream>
#include <string>
#include <vector>


ContextSV::ContextSV(Common common)
{
    this->common = common;
}

// Entry point
int ContextSV::run()
{
    // Call CNVs using the SNP positions
    std::cout << "Calling CNVs..." << std::endl;
    CNVCaller cnv_obj(this->common);
    CNVData cnv_calls = cnv_obj.run();

    // Call SVs from long read alignments and CNV calls
    // Return a map of SV type by start and end position
    // Key = [chromosome, SV start position], Value = [SV end position, SV type]
    std::cout << "Calling SVs..." << std::endl;
    SVCaller sv_obj(this->common);
    SVData sv_calls = sv_obj.run();

    // Classify SVs based on CNV calls
    std::cout << "Labeling CNVs..." << std::endl;
    this->labelCNVs(cnv_calls, sv_calls);

    // Write SV calls to file
    std::cout << "Writing SV calls to file..." << std::endl;
    std::string output_dir = this->common.getOutputDir();
    sv_calls.saveToVCF(output_dir);

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

        SVCandidate key = sv_call.first;
        
        // Get the SV coordinates
        std::string chr = std::get<0>(key);
        int start_pos = std::get<1>(key);
        int end_pos = std::get<2>(key);
        //int sv_type = std::get<3>(key);

        // Get CNV calls within the SV coordinate range and identify the most
        // common call
        int cnv_call = cnv_calls.getMostCommonCNV(chr, start_pos, end_pos);

        // Update the SV call's type
        sv_calls.updateSVType(chr, start_pos, end_pos, cnv_call);
    }
}