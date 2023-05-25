#include "common.h"

#include <stdio.h>
#include <string.h>

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

        // Get the start position
        else if (col == 1)
        {
            this->region_start = atoi(tok);
        }

        // Get the end position
        else if (col == 2)
        {
            this->region_end = atoi(tok);
        }

        tok = strtok(NULL, ":");
        col++;
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
