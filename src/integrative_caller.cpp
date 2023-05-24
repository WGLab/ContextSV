
#include "integrative_caller.h"
#include "snv_caller.h"
#include "cnv_caller.h"

#include <htslib/sam.h>
#include <iostream>
#include <string>


IntegrativeCaller::IntegrativeCaller()
= default;

/// Entry point
int IntegrativeCaller::run()
{
    // Check if the bam file uses chr prefix notation
    std::string filepath = bam_filepath;

    // Call SNVs
    //SNVCaller snv_obj;
    //snv_obj.run(filepath);

    // Call CNVs
    CNVCaller cnv_obj;
    cnv_obj.set_bam_filepath(filepath);
    cnv_obj.set_ref_filepath(ref_filepath);
    cnv_obj.set_output_dir(output_dir);
    cnv_obj.run(this->region_chr, this->region_start, this->region_end, this->window_size);

    return 0;
}

void IntegrativeCaller::set_bam_filepath(std::string bam_filepath)
{
    this->bam_filepath = bam_filepath;
}

void IntegrativeCaller::set_ref_filepath(std::string ref_filepath)
{
    this->ref_filepath = ref_filepath;
}

void IntegrativeCaller::set_output_dir(std::string output_dir)
{
    this->output_dir = output_dir;

    // Create the output directory
    std::string cmd = "mkdir -p " + output_dir;
    system(cmd.c_str());
}

void IntegrativeCaller::set_region(std::string region)
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
