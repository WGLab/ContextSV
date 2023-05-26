
#include "snv_caller.h"

#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <iostream>
#include <fstream>
# include <stdio.h>
#include <string>

SNVCaller::SNVCaller(Common common)
{
    this->common = common;
}

// Entry point
int SNVCaller::run()
{
    // Get file paths
    std::string bam_filepath = this->common.get_bam_filepath();
    std::string ref_filepath = this->common.get_ref_filepath();
    std::string output_dir = this->common.get_output_dir();

    // Use bcftools to call SNPs
    std::string cmd = "bcftools mpileup -Ou -f " + ref_filepath + " " + bam_filepath + " | bcftools call -mv -Ob -o " + output_dir + "/snps.bcf";

    // Run the command
    std::cout << "Calling SNPs..." << std::endl;
    std::cout << cmd << std::endl;
    int ret = system(cmd.c_str());
    if (ret != 0) {
        std::string err_str = "Error calling SNPs.";
        throw std::runtime_error(err_str);
    }

    // Filter the SNP BCF file
    cmd = "bcftools filter -i 'QUAL>20 && DP>10' " + output_dir + "/snps.bcf > " + output_dir + "/snps.filtered.bcf";
    std::cout << "Filtering SNP BCF file..." << std::endl;
    std::cout << cmd << std::endl;
    ret = system(cmd.c_str());
    if (ret != 0) {
        std::string err_str = "Error filtering SNP BCF file.";
        throw std::runtime_error(err_str);
    }

    // Read the SNP positions from the filtered BCF file into a vector
    std::cout << "Reading SNP positions from filtered BCF file..." << std::endl;
    htsFile* snp_bcf_file = bcf_open((output_dir + "/snps.filtered.bcf").c_str(), "r");
    if (snp_bcf_file == NULL) {
        std::string err_str = "Error opening SNP BCF file.";
        throw std::runtime_error(err_str);
    }

    // Read the header
    bcf_hdr_t* snp_bcf_header = bcf_hdr_read(snp_bcf_file);
    if (snp_bcf_header == NULL) {
        std::string err_str = "Error reading SNP BCF header.";
        throw std::runtime_error(err_str);
    }

    // Read the records
    bcf1_t* snp_bcf_record = bcf_init();
    if (snp_bcf_record == NULL) {
        std::string err_str = "Error initializing SNP BCF record.";
        throw std::runtime_error(err_str);
    }

    // Iterate over the records
    while (bcf_read(snp_bcf_file, snp_bcf_header, snp_bcf_record) == 0) {
        // Get the chromosome
        std::string chr = bcf_hdr_id2name(snp_bcf_header, snp_bcf_record->rid);

        // Get the position
        int pos = snp_bcf_record->pos + 1;

        // Print the position
        //std::cout << chr << ":" << pos << std::endl;

        // Add the position to the vector
        this->snp_positions.push_back(pos);
    }

    // Close the BCF file
    bcf_close(snp_bcf_file);

    // Save the SNP positions to CSV
    std::cout << "Saving SNP positions to CSV..." << std::endl;
    std::ofstream snp_csv_file;
    snp_csv_file.open(output_dir + "/snps.filtered.csv");
    for (int i = 0; i < this->snp_positions.size(); i++) {
        snp_csv_file << this->snp_positions[i] << std::endl;
    }
    snp_csv_file.close();

    return 0;
}

std::vector<int> SNVCaller::get_snp_positions()
{
    return this->snp_positions;
}
