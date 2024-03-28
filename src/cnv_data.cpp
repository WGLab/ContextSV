#include "cnv_data.h"

/// @cond
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>

#include "sv_types.h"
/// @endcond

// Include the SV types namespace
using namespace sv_types;

void CNVData::addCNVCall(std::string chr, int snp_pos, int cnv_type)
{
    // Add the CNV call to the map
    SNPLocation key(chr, snp_pos);
    this->cnv_calls[key] = cnv_type;
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
        snp_pos = std::stoi(pos_str);

        std::string skip_str;
        std::getline(iss, skip_str, '\t');
        std::getline(iss, skip_str, '\t');

        std::string cnv_type_str;
        std::getline(iss, cnv_type_str, '\t');
        cnv_type = std::stoi(cnv_type_str);

        // Add the CNV call to the map
        this->addCNVCall(chr, snp_pos, cnv_type);

        line_num++;
    }
    cnv_file.close();

    std::cout << "Loaded " << line_num << " CNV calls" << std::endl;
}
