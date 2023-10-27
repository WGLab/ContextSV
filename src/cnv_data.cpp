#include "cnv_data.h"

/// @cond
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
/// @endcond

void CNVData::addCNVCall(std::string chr, int snp_pos, int cnv_type)
{
    // Add the CNV call to the map
    SNPLocation key(chr, snp_pos);
    this->cnv_calls[key] = cnv_type;
}

int CNVData::getMostCommonCNV(std::string chr, int start, int end)
{
    // Get the majority CNV type within the SV region start and end positions
    // (0 = deletion, 1 = duplication, -1 = no CNV call)
    int dup_count = 0;
    int del_count = 0;
    int no_call_count = 0;
    int total_count = 0;
    //int sv_len = end - start + 1;
    
    //std::cout << "Checking for CNV calls in " << chr << ":" << start << "-" << end << " (SVLEN=" << sv_len << ")" << std::endl;

    for (int pos = start; pos <= end; pos++) {
        SNPLocation key(chr, pos);
        //std::cout << "Checking for CNV call at " << chr << ":" << pos <<
        //std::endl;
        
        if (this->cnv_calls.find(key) != this->cnv_calls.end()) {
            //std::cout << "State = " << this->cnv_calls[key] << std::endl;
            if (this->cnv_calls[key] == 5 || this->cnv_calls[key] == 6) {
                dup_count++;
            } else if (this->cnv_calls[key] == 1 || this->cnv_calls[key] == 2) {
                del_count++;
            } else {
                no_call_count++;
            }
            total_count++;
        }
    }

    // Check if the SV region is mostly covered by CNV calls (at least 50%) and
    // if the majority CNV type is an insertion or deletion
    int cnv_type = -1;
    if (total_count > 0) {
        if (dup_count > del_count && (double) dup_count / total_count > 0.5) {
            cnv_type = 1;
            //std::cout << "CNV type is DUP, SVLEN=" << sv_len << std::endl;
        } else if (del_count > dup_count && (double) del_count / total_count > 0.5) {
            cnv_type = 0;
            //std::cout << "CNV type is DEL, SVLEN=" << sv_len << std::endl;
        } else {
            //std::cout << "CNV type is no call, SVLEN=" << sv_len << std::endl;
        }
    }

    return cnv_type;
}

void CNVData::loadFromFile(std::string filepath)
{
    // Load CNV calls from file
    std::ifstream cnv_file(filepath);
    std::string line;
    std::string chr;
    int snp_pos;
    int cnv_type;

    // Check if the file was opened successfully
    if (!cnv_file.is_open()) {
        std::cerr << "Error: Could not open CNV file " << filepath << std::endl;
        exit(1);
    }

    // Skip the first line (header)
    std::getline(cnv_file, line);

    // Read the file line by line
    int line_num = 1;
    while (std::getline(cnv_file, line)) {

        // Parse the line
        std::istringstream iss(line);

        // Get columns 1, 2, and 5 (chr, pos, cnv_type)
        std::string chr;
        std::getline(iss, chr, '\t');

        std::string pos_str;
        std::getline(iss, pos_str, '\t');
        int snp_pos = std::stoi(pos_str);

        std::string skip_str;
        std::getline(iss, skip_str, '\t');
        std::getline(iss, skip_str, '\t');

        std::string cnv_type_str;
        std::getline(iss, cnv_type_str, '\t');
        int cnv_type = std::stoi(cnv_type_str);

        // Add the CNV call to the map
        this->addCNVCall(chr, snp_pos, cnv_type);

        line_num++;
    }
    cnv_file.close();

    std::cout << "Loaded " << line_num << " CNV calls" << std::endl;
}
