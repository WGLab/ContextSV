#include "contextsv.h"
#include "cnv_caller.h"
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
    // Get the reference genome
    FASTAQuery ref_genome = this->input_data->getRefGenome();

    // Call SVs from long read alignments:
    std::cout << "Running alignment-based SV calling..." << std::endl;
    auto start_sv = std::chrono::high_resolution_clock::now();
    SVCaller sv_caller(*this->input_data);
    SVData sv_calls = sv_caller.run();
    auto end_sv = std::chrono::high_resolution_clock::now();

    // Format and print the time taken to call SVs
    std::string elapsed_time = getElapsedTime(start_sv, end_sv);
    std::cout << "Alignment-based SV calling complete. Found " << sv_calls.size() << " total SVs. Time taken (h:m:s) = " << elapsed_time << std::endl;

    // Classify SVs based on SNP CNV predictions if enabled
    if (this->input_data->getDisableSNPCNV() == false) {

        std::cout << "Running SNP CNV calling..." << std::endl;

        // Check if a file with CNV data was provided
        CNVData cnv_calls;
        if (this->input_data->getCNVFilepath() != "") {
            // Load CNV data
            std::cout << "Loading CNV data..." << std::endl;
            cnv_calls.loadFromFile(this->input_data->getCNVFilepath());
        } else {
            // Call CNVs at SNP positions
            auto start_cnv = std::chrono::high_resolution_clock::now();
            CNVCaller cnv_caller(*this->input_data);
            cnv_caller.run(sv_calls);
            auto end_cnv = std::chrono::high_resolution_clock::now();

            // Format and print the time taken to call CNVs
            elapsed_time = getElapsedTime(start_cnv, end_cnv);
            std::cout << "SNP CNV calling complete. Time taken (h:m:s) = " << elapsed_time << std::endl;
        }

        std::cout << "Running SV CNV labeling from SNP predictions..." << std::endl;
        auto start_label = std::chrono::high_resolution_clock::now();
        //this->labelCNVs(cnv_calls, sv_calls);
        auto end_label = std::chrono::high_resolution_clock::now();
        elapsed_time = getElapsedTime(start_label, end_label);
        std::cout << "SV CNV labeling complete. Time taken (h:m:s) = " << elapsed_time << std::endl;
    }

    // Write SV calls to file
    std::string output_dir = this->input_data->getOutputDir();
    std::cout << "Writing SV calls to file " << output_dir << "/sv_calls.vcf..." << std::endl;
    auto start_write = std::chrono::high_resolution_clock::now();
    sv_calls.saveToVCF(ref_genome, output_dir);
    auto end_write = std::chrono::high_resolution_clock::now();
    elapsed_time = getElapsedTime(start_write, end_write);
    std::cout << "SV calls written to file. Time taken (h:m:s) = " << elapsed_time << std::endl;

    return 0;
}

// Label SVs based on CNV calls
void ContextSV::labelCNVs(CNVData& cnv_calls, SVData& sv_calls)
{
    // Iterate over SV calls
    for (auto const& sv_call : sv_calls) {

        SVCandidate candidate = sv_call.first;
        
        // Get the SV coordinates
        std::string chr = std::get<0>(candidate);
        int start_pos = std::get<1>(candidate);
        int end_pos = std::get<2>(candidate);

        // Get CNV calls within the SV coordinate range and identify the most
        // common call
        std::tuple<int, std::string> cnv_info = cnv_calls.getMostCommonCNV(chr, start_pos, end_pos);

        // Get the CNV type and genotype
        int cnv_type = std::get<0>(cnv_info);
        std::string cnv_genotype = std::get<1>(cnv_info);

        // Update the SV call's type if the CNV call is not unknown
        if (cnv_type != UNKNOWN) {
            sv_calls.updateSVType(candidate, cnv_type, "SNPCNV");
        }

        // Update the SV call's genotype
        sv_calls.updateGenotype(candidate, cnv_genotype);
    }
}
