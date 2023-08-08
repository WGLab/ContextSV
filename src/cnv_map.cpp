#include "cnv_map.h"

#include <iostream>

void CNVMap::addCNVCall(std::string chr, int snp_pos, int cnv_type)
{
    std::pair<std::string, int> key(chr, snp_pos);
    this->cnv_calls[key] = cnv_type;
    //std::cout << "Added SNP CNV call at " << chr << ":" << snp_pos << std::endl;
}

std::map<std::pair<std::string, int>, int> CNVMap::getCNVCalls()
{
    return this->cnv_calls;
}

int CNVMap::getSVType(std::string chr, int start, int end)
{
    // Get the majority CNV type within the SV region start and end positions
    // (1=INS, 2=DEL)
    int ins_count = 0;
    int del_count = 0;
    int no_call_count = 0;
    int total_count = 0;
    
    //std::cout << "Checking for CNV calls in " << chr << ":" << start << "-" << end << std::endl;

    for (int pos = start; pos <= end; pos++) {
        std::pair<std::string, int> key(chr, pos);
        //std::cout << "Checking for CNV call at " << chr << ":" << pos <<
        //std::endl;
        
        if (this->cnv_calls.find(key) != this->cnv_calls.end()) {
            //std::cout << "Found CNV call at " << chr << ":" << pos << std::endl;
            if (this->cnv_calls[key] == 5 || this->cnv_calls[key] == 6) {
                ins_count++;
            } else if (this->cnv_calls[key] == 1 || this->cnv_calls[key] == 2) {
                del_count++;
            } else {
                no_call_count++;
            }
            total_count++;
        }
    }

    // Summarize the CNV calls
    // std::cout << "Total CNV calls: " << total_count << std::endl;
    // std::cout << "INS calls: " << ins_count << std::endl;
    // std::cout << "DEL calls: " << del_count << std::endl;
    // std::cout << "No call: " << no_call_count << std::endl;


    // Check if the SV region is covered by CNV calls
    if ((total_count > 0) && (ins_count > no_call_count || del_count > no_call_count)) {

        // Check if the majority CNV type is an insertion
        if (ins_count > del_count) {
            return 1;
        }

        // Check if the majority CNV type is a deletion
        if (del_count > ins_count) {
            return 2;
        }
    }

    // // Print the percentage of CNV calls
    // std::cout << "INS calls: " << ins_count << " (" << (ins_count / total_count) * 100 << "%)" << std::endl;
    // std::cout << "DEL calls: " << del_count << " (" << (del_count / total_count) * 100 << "%)" << std::endl;
    // std::cout << "No call: " << no_call_count << " (" << (no_call_count / total_count) * 100 << "%)" << std::endl;

    // Return 0 if the SV region is not covered by CNV calls
    return 0;
}
