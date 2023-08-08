#ifndef CNV_MAP_H
#define CNV_MAP_H

#include <string>
#include <map>

class CNVMap {
    private:
        std::map<std::pair<std::string, int>, int> cnv_calls;

    public:
        void addCNVCall(std::string chr, int snp_pos, int cnv_type);
        std::map<std::pair<std::string, int>, int> getCNVCalls();
        int getSVType(std::string chr, int start, int end);
};

#endif // CNV_MAP_H
