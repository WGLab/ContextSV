
#include "integrative_caller.h"
#include "snv_caller.h"
#include "cnv_caller.h"
#include "common.h"

#include <htslib/sam.h>
#include <iostream>
#include <string>


IntegrativeCaller::IntegrativeCaller(Common common)
{
    this->common = common;
}

/// Entry point
int IntegrativeCaller::run()
{
    // Call SNVs
    // SNVCaller snv_obj(this->common);
    // snv_obj.run();

    // Get the positions of SNPs in the region
    //snv_obj.run(filepath);

    // Print the size of the SNP positions vector
    // std::cout << "SNP positions vector size = " << snv_obj.get_snp_positions().size() << std::endl;

    // Call CNVs
    //CNVCaller cnv_obj(this->common, snv_obj.get_snp_positions());
    std::vector<int> snp_positions;
    CNVCaller cnv_obj(this->common, snp_positions);
    cnv_obj.run();

    return 0;
}
