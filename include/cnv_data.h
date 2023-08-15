#ifndef CNV_DATA_H
#define CNV_DATA_H

#include <string>
#include <map>

// Type aliases
class CNVData;
        // SNP location (chr, snp_pos)
        typedef std::pair<std::string, int> SNPLocation;

        // CNV map (SNP location, CNV type)
        typedef std::map<SNPLocation, int> SNPToCNVMap;

class CNVData {
    private:
        SNPToCNVMap cnv_calls;

    public:
        void addCNVCall(std::string chr, int snp_pos, int cnv_type);
        int getMostCommonCNV(std::string chr, int start, int end);
};

#endif // CNV_DATA_H
