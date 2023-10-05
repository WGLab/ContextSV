#include "cnv_data.h"

/// @cond
#include <iostream>
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
    int sv_len = end - start + 1;
    
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
