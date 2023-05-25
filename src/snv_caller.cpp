
#include "snv_caller.h"

#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <iostream>
# include <stdio.h>
#include <string>

SNVCaller::SNVCaller(Common common)
{
    this->common = common;
}

// Entry point
int SNVCaller::run()
{
	// Get the SNP VCF file path
    std::string snp_vcf_filepath = this->common.get_snp_vcf_filepath();

    // Read the SNP VCF file
    htsFile *snp_vcf_file = bcf_open(snp_vcf_filepath.c_str(), "r");
    bcf_hdr_t *snp_vcf_header = bcf_hdr_read(snp_vcf_file);

    // Get the positions of all SNPs and add them to a vector
    std::vector<int> snp_positions;
    bcf1_t *snp_vcf_record = bcf_init();
    while (bcf_read(snp_vcf_file, snp_vcf_header, snp_vcf_record) == 0) {
        int snp_position = snp_vcf_record->pos + 1;
        snp_positions.push_back(snp_position);

        // Get the B-allele count
        //int b_allele_count = bcf_get_format_int32(snp_vcf_header, snp_vcf_record, "BAC");

    }
    // Close the SNP VCF file
    bcf_close(snp_vcf_file);

    // Update the SNP positions
    this->snp_positions = snp_positions;

    return 0;
}

std::vector<int> SNVCaller::get_snp_positions()
{
    return this->snp_positions;
}
