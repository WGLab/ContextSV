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

        // Define a map of CNV genotypes by HMM predicted state.
        // Each of the 6 state predictions corresponds to a copy number state:
        // 1: 0/0 (Two copy loss: homozygous deletion, GT: 0/0)
        // 2: 1/0 (One copy loss: heterozygous deletion, GT: 0/1)
        // 3: 1/1 (Normal diploid: no copy number change, GT: 1/1)
        // 4: 1/1 (Copy neutral LOH: no copy number change, GT: 1/1)
        // 5: 2/1 (One copy gain: heterozygous duplication, GT: 1/2)
        // 6: 2/2 (Two copy gain: homozygous duplication, GT: 2/2)
        std ::map<int, std::string> cnv_genotype_map = {
            {1, "0/0"},
            {2, "0/1"},
            {3, "1/1"},
            {4, "1/1"},
            {5, "1/2"},
            {6, "2/2"}
        };

    public:
        // Add a CNV call to the map
        void addCNVCall(std::string chr, int snp_pos, int cnv_type);

        // Get the most common CNV type within the SV region start and end positions
        std::tuple<int, std::string> getMostCommonCNV(std::string chr, int start, int end);

        // Load CNV calls from file
        void loadFromFile(std::string filepath);
};

#endif // CNV_DATA_H
