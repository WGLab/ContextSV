#include "input_data.h"

#include <stdio.h>
#include <string.h>
#include <iostream>

#define BUFFER_SIZE 1024

// Constructor
InputData::InputData()
{
    // this->region_set = false;
    // this->window_size = 1000;
}

std::string InputData::getBAMFilepath()
{
    return this->bam_filepath;
}

void InputData::setBAMFilepath(std::string filepath)
{
    this->bam_filepath = filepath;
}

void InputData::setRefGenome(std::string fasta_filepath)
{
    // Set the reference genome
    this->fasta_query.setFilepath(fasta_filepath);
}

FASTAQuery InputData::getRefGenome()
{
    return this->fasta_query;
}

std::string InputData::getOutputDir()
{
    return this->output_dir;
}

void InputData::setOutputDir(std::string dirpath)
{
    this->output_dir = dirpath;

    // Create the output directory
    std::string cmd = "mkdir -p " + output_dir;
    system(cmd.c_str());
}

std::string InputData::getRegion()
{
    return this->region;
}

void InputData::setRegion(std::string region)
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

int InputData::getWindowSize()
{
    return this->window_size;
}

void InputData::setWindowSize(int window_size)
{
    this->window_size = window_size;
}

std::string InputData::getSNPFilepath()
{
    return this->snp_vcf_filepath;
}

void InputData::setSNPFilepath(std::string filepath)
{
    this->snp_vcf_filepath = filepath;
}

std::string InputData::getRegionChr()
{
    return this->region_chr;
}

int InputData::getRegionStart()
{
    return this->region_start;
}

int InputData::getRegionEnd()
{
    return this->region_end;
}

bool InputData::getRegionSet()
{
    return this->region_set;
}
