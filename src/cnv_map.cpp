#include "cnv_map.h"



void CNVMap::addCNVCall(char * chr, int snp_pos, int cnv_type)
{
    std::pair<char *, int> key(chr, snp_pos);
    this->cnv_calls[key] = cnv_type;
}

std::map<std::pair<char *, int>, int> CNVMap::getCNVCalls()
{
    return this->cnv_calls;
}
