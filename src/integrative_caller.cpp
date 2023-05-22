
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
    bool uses_chr_prefix;
    int exit_code = bamHasChrPrefix(filepath, uses_chr_prefix);

    // Call SNVs
    //SNVCaller snv_obj;
    //snv_obj.run(filepath);

    // Call CNVs
    CNVCaller cnv_obj;
    cnv_obj.set_bam_filepath(filepath);
    cnv_obj.set_ref_filepath(ref_filepath);
    cnv_obj.set_output_dir(output_dir);
    cnv_obj.run(this->region_chr, this->region_start, this->region_end);

    return 0;
}

/// Check if the bam file uses chr prefix notation
int IntegrativeCaller::bamHasChrPrefix(std::string filepath, bool& uses_chr_prefix)
{
    samFile* bam_file = sam_open(filepath.c_str(), "r");
    if (bam_file == nullptr) {
        std::cerr << "Error: failed to open " << bam_file << "\n";
        return 1;
    }

    bam_hdr_t* header = sam_hdr_read(bam_file);
    if (header == nullptr) {
        std::cerr << "Error: failed to read header from " << filepath << "\n";
        sam_close(bam_file);
        return 1;
    }

    uses_chr_prefix = false;
    for (int i = 0; i < header->n_targets; i++) {
        const char* name = header->target_name[i];
        if (name[0] == 'c' && name[1] == 'h' && name[2] == 'r') {
            uses_chr_prefix = true;
            break;
        }
    }

    std::cout << "The bam file " << filepath << " uses " << (uses_chr_prefix ? "chr" : "") << " notation\n";

    bam_hdr_destroy(header);
    sam_close(bam_file);
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
