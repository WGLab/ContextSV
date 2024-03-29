#ifndef CNV_DATA_H
#define CNV_DATA_H

#include "types.h"

/// @cond
#include <string>
/// @endcond


class CNVData {
    private:
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
