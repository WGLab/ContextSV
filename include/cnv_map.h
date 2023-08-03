#ifndef CNV_MAP_H
#define CNV_MAP_H

#include <map>

class CNVMap {
    private:
        std::map<std::pair<char *, int>, int> cnv_calls;

    public:
        void addCNVCall(char * chr, int snp_pos, int cnv_type);
        std::map<std::pair<char *, int>, int> getCNVCalls();
};

#endif // CNV_MAP_H
