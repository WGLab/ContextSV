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

        // Get the start position
        else if (col == 1)
        {
            this->region_start = atoi(tok);
            std::cout << "Region start = " << this->region_start << std::endl;
        }

        // Get the end position
        else if (col == 2)
        {
            this->region_end = atoi(tok);
            std::cout << "Region end = "  << this->region_end << std::endl;
        }

        tok = strtok(NULL, ":");
        col++;
    }

    // Calculate the chromosome length
    std::cout << "Getting chromosome length for " << this->region_chr << std::endl;
    std::string input_filepath = this->get_bam_filepath();
    std::string target_chr = this->get_region_chr();
    std::string chr_length_cmd = "samtools view -H " + input_filepath + " | grep \"@SQ\" | grep " + target_chr + " | cut -f 3 -d ':' | cut -f 2 -d '-'";
    FILE *chr_length_pipe = popen(chr_length_cmd.c_str(), "r");
    char chr_length_buffer[BUFFER_SIZE];
    fgets(chr_length_buffer, BUFFER_SIZE, chr_length_pipe);
    pclose(chr_length_pipe);
    int chr_length = atoi(chr_length_buffer);
    std::cout << "Chromosome length = " << chr_length << std::endl;
    this->chr_length = chr_length;
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

int Common::get_chr_length()
{
    return this->chr_length;
}
