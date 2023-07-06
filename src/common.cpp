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

void Common::print_progress(int progress, int total)
{
    // Get the percentage
    float percent = (float)progress / (float)total * 100.0;

    // Get the number of hashes
    int num_hashes = (int)(percent / 2.0);

    // Print the progress bar
    printf("\r[");
    for (int i = 0; i < num_hashes; i++)
    {
        printf("#");
    }
    for (int i = 0; i < 50 - num_hashes; i++)
    {
        printf(" ");
    }
    printf("] %3.2f%%", percent);
    fflush(stdout);

    // Print a new line if finished
    if (progress == total)
    {
        printf("\n");
    }
}
