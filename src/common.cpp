#include "common.h"

#include <stdio.h>
#include <string.h>
#include <iostream>

#define BUFFER_SIZE 1024

std::string Common::get_bam_filepath()
{
    return this->bam_filepath;
}

void Common::set_bam_filepath(std::string filepath)
{
    this->bam_filepath = filepath;
}

std::string Common::get_ref_filepath()
{
    return this->ref_filepath;
}

void Common::set_ref_filepath(std::string filepath)
{
    this->ref_filepath = filepath;
}

std::string Common::get_output_dir()
{
    return this->output_dir;
}

void Common::set_output_dir(std::string dirpath)
{
    this->output_dir = dirpath;

    // Create the output directory
    std::string cmd = "mkdir -p " + output_dir;
    system(cmd.c_str());
}

std::string Common::get_region()
{
    return this->region;
}

void Common::set_region(std::string region)
{
    this->region = region;

    // Parse the region
    char *tok = strtok((char *)region.c_str(), ":");
    int col = 0;
    while (tok != NULL)
    {
        // Get the chromosome
        if (col == 0)
        {
            this->region_chr = tok;
        }

        // Get the start and end positions
        else if (col == 1)
        {
            // Check if empty
            if (strcmp(tok, "") == 0)
            {
                break;
            }

            // Split the start and end positions
            char *start_tok = strtok(tok, "-");
            char *end_tok = strtok(NULL, "-");

            // Get the start position
            if (start_tok != NULL)
            {
                this->region_start = atoi(start_tok);
            }

            // Get the end position
            if (end_tok != NULL)
            {
                this->region_end = atoi(end_tok);
            }
        }
        tok = strtok(NULL, ":");
        col++;
    }

    // // Calculate the chromosome length
    // std::cout << "Getting chromosome length for " << this->region_chr << std::endl;
    // std::string input_filepath = this->get_bam_filepath();
    // std::string target_chr = this->get_region_chr();
    // std::string chr_length_cmd = "samtools view -H " + input_filepath + " | grep \"@SQ\" | grep " + target_chr + " | cut -f 3 -d ':' | cut -f 2 -d '-'";
    // FILE *chr_length_pipe = popen(chr_length_cmd.c_str(), "r");
    // char chr_length_buffer[BUFFER_SIZE];
    // fgets(chr_length_buffer, BUFFER_SIZE, chr_length_pipe);
    // pclose(chr_length_pipe);
    // int chr_length = atoi(chr_length_buffer);

    // // Print the chromosome length
    // std::cout << "Chromosome length = " << chr_length << std::endl;
    // this->chr_length = chr_length;

    // Print the region
    if (this->region_start == 0 && this->region_end == 0)
    {
        this->region_set = false;
        std::cout << "Parsed region = " << this->region_chr << std::endl;
    } else {
        this->region_set = true;
        std::cout << "Parsed region = " << this->region_chr << ":" << this->region_start << "-" << this->region_end << std::endl;
    }
}

int Common::get_window_size()
{
    return this->window_size;
}

void Common::set_window_size(int window_size)
{
    this->window_size = window_size;
}

std::string Common::get_snp_vcf_filepath()
{
    return this->snp_vcf_filepath;
}

void Common::set_snp_vcf_filepath(std::string filepath)
{
    this->snp_vcf_filepath = filepath;

    // // Create a VCF filepath of filtered SNPs
    // std::string filtered_snp_vcf_filepath = this->output_dir + "/filtered_snps.vcf";

    // std::cout << "Parsing SNPs from " << this->snp_vcf_filepath << std::endl;

    // // Filter variants by depth and quality and SNPs only
    // std::string cmd = "bcftools view -r " + this->region + " -v snps -i 'QUAL > 30 && DP > 10 && FILTER = \"PASS\"' " + this->snp_vcf_filepath + " > " + filtered_snp_vcf_filepath;
    // std::cout << "Command: " << cmd << std::endl;
    // system(cmd.c_str());

    // std::cout << "Filtered SNPs written to " << filtered_snp_vcf_filepath << std::endl;

    // // Extract all BAFs from the filtered SNPs
    // std::string baf_filepath = this->output_dir + "/filtered_bafs.csv";
    // cmd = "bcftools query -f '%CHROM,%POS,[%VAF]\n' " + filtered_snp_vcf_filepath + " > " + baf_filepath;
    // std::cout << "Command: " << cmd << std::endl;
    // system(cmd.c_str());

    // std::cout << "BAFs written to " << baf_filepath << std::endl;
}

std::string Common::get_region_chr()
{
    return this->region_chr;
}

int Common::get_region_start()
{
    return this->region_start;
}

int Common::get_region_end()
{
    return this->region_end;
}

bool Common::get_region_set()
{
    return this->region_set;
}
