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
        // Define constants for CNV types
        static const int DEL = 0;  // Deletion
        static const int DUP = 1;  // Duplication
        static const int NEUT = 2;  // Copy neutral
        static const int UNKNOWN = -1;  // Unknown
        SNPToCNVMap cnv_calls;  // Map of SNP positions to CNV types

    public:
        // Add a CNV call to the map
        void addCNVCall(std::string chr, int snp_pos, int cnv_type);

        // Get the most common CNV type within the SV region start and end positions
        int getMostCommonCNV(std::string chr, int start, int end);

        // Load CNV calls from file
        void loadFromFile(std::string filepath);
};

#endif // CNV_DATA_H
