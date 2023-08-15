#include "sv_data.h"

#include <iostream>
#include <htslib/vcf.h>  // For saving SV calls to VCF
#include <htslib/vcfutils.h>  // For saving SV calls to VCF
#include <htslib/sam.h>  // For reading BAM files

void SVData::addSVCall(std::string chr, int start, int end, int sv_type)
{

    // Add the SV call to the map of candidate locations
    SVCandidate key(chr, start, end);
    if (this->sv_calls.find(key) != this->sv_calls.end()) {
        this->sv_calls[key].first = sv_type;
        this->sv_calls[key].second++;
    } else {
        this->sv_calls[key] = SVInfo(sv_type, 1);
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
            this->sv_calls[key].first = value.first;
            this->sv_calls[key].second += value.second;
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
        this->sv_calls[key].first = sv_type;
    }
}

void SVData::saveToVCF(std::string filename)
{
    // Open the VCF file for writing
    htsFile* vcf_file = bcf_open(filename.c_str(), "w");
    if (vcf_file == NULL) {
        std::cerr << "Error: could not open VCF file " << filename << "\n";
        exit(1);
    }

    // Create the VCF header
    bcf_hdr_t* vcf_header = bcf_hdr_init("w");
    if (vcf_header == NULL) {
        std::cerr << "Error: could not create VCF header\n";
        exit(1);
    }

    // Add the contig information to the VCF header (TODO: add all contigs)
    // std::string contig = "##contig=<ID=" + chr + ">";
    // bcf_hdr_append(vcf_header, contig.c_str());

    // Add the INFO field to the VCF header
    std::string info = "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">";
    bcf_hdr_append(vcf_header, info.c_str());
    info = "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of structural variant\">";
    bcf_hdr_append(vcf_header, info.c_str());
    info = "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of structural variant\">";
    bcf_hdr_append(vcf_header, info.c_str());
    info = "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Read depth at this position\">";
    bcf_hdr_append(vcf_header, info.c_str());

    // Add the sample to the VCF header
    std::string sample = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE";
    bcf_hdr_append(vcf_header, sample.c_str());

    // Add the SV calls to the VCF file
    for (auto const& sv_call : this->sv_calls) {

        std::cout << "Saving SV call to VCF file\n";

        // Get the SV candidate location
        SVCandidate key = sv_call.first;
        std::string chr = std::get<0>(key);
        int start = std::get<1>(key);
        int end = std::get<2>(key);

        // Get the SV type and read depth
        SVInfo value = sv_call.second;
        int sv_type = value.first;
        int read_depth = value.second;

        // Create a VCF header line
        bcf_hdr_t* vcf_header_line = bcf_hdr_init("w");
        if (vcf_header_line == NULL) {
            std::cerr << "Error: could not create VCF header line\n";
            exit(1);
        }

        // Write the VCF header line to the VCF file
        if (bcf_hdr_write(vcf_file, vcf_header_line) != 0) {
            std::cerr << "Error: could not write VCF header line\n";
            exit(1);
        }

        // Create the VCF record
        bcf1_t* vcf_record = bcf_init1();
        if (vcf_record == NULL) {
            std::cerr << "Error: could not create VCF record\n";
            exit(1);
        }

        std::cout << "Adding sample information to VCF record\n";

        // Add the CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, and FORMAT
        // fields to the VCF record
        // bcf_update_info_string(vcf_header, vcf_record, "CHROM", chr.c_str());
        // bcf_update_info_int32(vcf_header, vcf_record, "POS", &start, 1);
        // bcf_update_info_string(vcf_header, vcf_record, "ID", ".");
        // bcf_update_info_string(vcf_header, vcf_record, "REF", ".");
        // bcf_update_info_string(vcf_header, vcf_record, "ALT", ".");
        // bcf_update_info_string(vcf_header, vcf_record, "QUAL", ".");
        // bcf_update_info_string(vcf_header, vcf_record, "FILTER", ".");
        // bcf_update_info_string(vcf_header, vcf_record, "INFO", ".");
        // bcf_update_info_string(vcf_header, vcf_record, "FORMAT", ".");

        std::cout << "Adding other info to VCF record\n";

        // // Set the end position
        // std::cout << "End position: " << end << std::endl;
        // int end_test = 1000;
        // int result = bcf_update_info_int32(vcf_header, vcf_record, "END", &end_test, 1);
        // if (result != 0) {
        //     std::cerr << "Error: could not set end position\n";
        //     exit(1);
        // }
        // //bcf_update_info_int32(vcf_header, vcf_record, "END", &end, 1);
        // std::cout << "End position set" << std::endl;

        // Set the SV type
        // Create a string to hold the SV type
        std::string sv_type_string;
        if (sv_type == 0) {
            sv_type_string = "DEL";  // Deletion
        } else if (sv_type == 1) {
            sv_type_string = "DUP";  // Duplication
        } else if (sv_type == 2) {
            sv_type_string = "INV";  // Inversion
        } else if (sv_type == 3) {
            sv_type_string = "INS";  // Insertion
        } else if (sv_type == 4) {
            sv_type_string = "BND";  // Translocation
        } else {
            sv_type_string = "NA";
        }

        // // Get the INFO SVTYPE field
        // bcf_info_t* info_field = bcf_get_info(vcf_header, "INFO/SVTYPE", NULL, NULL);


        // // Create a string array to hold the SV type
        // char* sv_type_array = new char[sv_type_string.length() + 1];

        // // Update the SVTYPE field
        // info_field->vptr = sv_type_array;
        
        // Set the read depth
        // Create an integer array to hold the read depth
        int* read_depth_array = new int[1];
        read_depth_array[0] = read_depth;
        std::cout << "[test] Read depth: " << read_depth << "\n";
        int result = bcf_update_info_int32(vcf_header, vcf_record, "INFO/DP", read_depth_array, 1);
        std::cout << "Result: " << result << "\n";

        // Clean up
        std::cout << "Cleaning up\n";
        delete[] read_depth_array;
        std::cout << "Cleaned up\n";
        if (result != 0) {
            std::cerr << "Error: could not set read depth\n";
            exit(1);
        }
        std::cout << "Read depth set\n";

        // Set the SV length
        int sv_length = end - start;
        bcf_update_info_int32(vcf_header, vcf_record, "SVLEN", &sv_length, 1);

        // Set the SV type
        if (sv_type == 1) {
            bcf_update_info_string(vcf_header, vcf_record, "SVTYPE", "INS");
        } else if (sv_type == 2) {
            bcf_update_info_string(vcf_header, vcf_record, "SVTYPE", "DEL");
        } else {
            bcf_update_info_string(vcf_header, vcf_record, "SVTYPE", "UNK");
        }

        // Write the VCF record to the VCF file
        if (bcf_write(vcf_file, vcf_header, vcf_record) != 0) {
            std::cerr << "Error: could not write VCF record\n";
            exit(1);
        }

        // Free the VCF record
        bcf_destroy(vcf_record);

        std::cout << "Wrote VCF record for " << chr << ":" << start << "-" << end << "\n";
    }

    // Close the VCF file
    bcf_close(vcf_file);

    // Free the VCF header
    bcf_hdr_destroy(vcf_header);

    // Free the VCF file
    hts_close(vcf_file);
}
