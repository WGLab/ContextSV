#ifndef CNV_DATA_H
#define CNV_DATA_H

#include "types.h"

#include <string>
#include <map>
#include <utility>  // std::pair


class CNVData {
    private:
        SNPToCNVMap cnv_calls;

    public:
        void addCNVCall(std::string chr, int snp_pos, int cnv_type);
        int getMostCommonCNV(std::string chr, int start, int end);
};

#endif // CNV_DATA_H
