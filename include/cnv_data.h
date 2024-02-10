#ifndef CNV_DATA_H
#define CNV_DATA_H

/// @cond
#include <string>
#include <set>
#include <map>
/// @endcond

// CNV candidate location map
// (chr, snp_pos) : cnv_type

using SNPLocation = std::pair<std::string, int64_t>;
using SNPToCNVMap = std::map<SNPLocation,  int>;

class CNVData {
    private:
        SNPToCNVMap cnv_calls;  // Map of SNP positions to CNV types

    public:
        // Add a CNV call to the map
        void addCNVCall(std::string chr, int snp_pos, int cnv_type);

        // Get the most common CNV type within the SV region start and end positions
        std::tuple<int, std::string> getMostCommonCNV(std::string chr, int start, int end);

        // Load CNV calls from file
        void loadFromFile(std::string filepath);
};

#endif // CNV_DATA_H
