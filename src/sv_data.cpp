#include "sv_data.h"

#include <iostream>
#include <fstream>

void SVData::addSVCall(std::string chr, int start, int end, int sv_type)
{

    // Add the SV call to the map of candidate locations
    SVCandidate key(chr, start, end);
    if (this->sv_calls.find(key) != this->sv_calls.end()) {
        std::get<0>(this->sv_calls[key]) = sv_type;
        std::get<1>(this->sv_calls[key]) += 1;
        std::get<2>(this->sv_calls[key]) = "";
        std::get<3>(this->sv_calls[key]) = "";
    } else {
        this->sv_calls[key] = SVInfo(sv_type, 1, "", "");
    }
}

void SVData::addSVCalls(SVData sv_calls)
{
    // Iterate over the SV calls
    for (auto const& sv_call : sv_calls) {

        // Add the SV call to the map
        SVCandidate key = sv_call.first;
        SVInfo value = sv_call.second;
        if (this->sv_calls.find(key) != this->sv_calls.end()) {
            std::get<0>(this->sv_calls[key]) = std::get<0>(value);
            std::get<1>(this->sv_calls[key]) += std::get<1>(value);
            std::get<2>(this->sv_calls[key]) = "";
            std::get<3>(this->sv_calls[key]) = "";
        } else {
            this->sv_calls[key] = value;
        }
    }
}

void SVData::updateSVType(std::string chr, int start, int end, int sv_type)
{
    // Update the SV type for a given SV candidate
    SVCandidate key(chr, start, end);
    if (this->sv_calls.find(key) != this->sv_calls.end()) {
        std::get<0>(this->sv_calls[key]) = sv_type;
    }
}

void SVData::saveToVCF(std::string output_dir)
{

    // Create a VCF file for the SV calls
    std::string output_vcf = output_dir + "/sv_calls.vcf";
    std::ofstream output_stream(output_vcf);

    // Write the header
    output_stream << "##fileformat=VCFv4.2" << std::endl;
    output_stream << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">" << std::endl;
    output_stream << "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">" << std::endl;
    output_stream << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth at the locus\">" << std::endl;
    output_stream << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" << std::endl;

    // Iterate over the SV calls
    std::cout << "Saving SV calls to " << output_vcf << "..." << std::endl;
    for (auto const& sv_call : this->sv_calls) {

        // Get the SV candidate and SV info
        SVCandidate key = sv_call.first;
        SVInfo value = sv_call.second;

        // Get the CHROM, POS, and END
        std::string chr = std::get<0>(key);
        int pos = std::get<1>(key);
        int end = std::get<2>(key);

        // Get the SV type and depth
        int sv_type = std::get<0>(value);
        int depth = std::get<1>(value);

        // Get the SV length
        int sv_len = end - pos;

        // Write the SV call to the file
        output_stream << chr << "\t" << pos << "\t" << "." << "\t" << "." << "\t" << "<" << sv_type << ">" << "\t" << "." << "\t" << "." << "\t" << "SVTYPE=" << sv_type << ";SVLEN=" << sv_len << ";DP=" << depth << std::endl;
    }

    // Close the output stream
    output_stream.close();
    std::cout << "Saved SV calls to " << output_vcf << std::endl;
}
