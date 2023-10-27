#include "input_data.h"

/// @cond
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <sstream>
/// @endcond

#define BUFFER_SIZE 1024

// Constructor
InputData::InputData()
{
    this->bam_filepath = "";
    this->ref_filepath = "";
    this->snp_vcf_filepath = "";
    this->pfb_filepath = "";
    this->output_dir = "";
    this->region = "";
    this->window_size = 10000;
    this->region_chr = "";
    this->region_start = 0;
    this->region_end = 0;
    this->region_set = false;
    this->thread_count = 1;
    this->hmm_filepath = "data/wgs.hmm";
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

    // Check if a valid chromosome was parsed
    if (this->region_chr == "")
    {
        std::cerr << "Error: Region chromosome not set" << std::endl;
        exit(1);
    }

    // Check if only a chromosome was parsed
    if (this->region_start == 0 && this->region_end == 0)
    {
        // Use the entire chromosome as the region
        this->region_set = false;
        std::cout << "Parsed region = " << this->region_chr << std::endl;
    } else {
        // Check if a valid chromosome start and end position were parsed
        if (this->region_start == 0 || this->region_end == 0 || this->region_start > this->region_end)
        {
            std::cerr << "Error: Region start and end positions not set" << std::endl;
            exit(1);
        } else {
            // Set the region
            this->region_set = true;
            std::cout << "Parsed region = " << this->region_chr << ":" << this->region_start << "-" << this->region_end << std::endl;
        }
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

void InputData::setChrCov(std::string chr_cov)
{
    // Update the chromosome coverage map if the string is not empty
    if (chr_cov != "")
    {
        // Split the string by commas
        std::istringstream ss(chr_cov);
        std::string token;
        std::vector<std::string> chr_cov_pairs;

        while (std::getline(ss, token, ','))
        {
            chr_cov_pairs.push_back(token);
        }

        // Iterate over the pairs
        for (auto const &pair : chr_cov_pairs)
        {
            // Split the pair by colon
            std::istringstream ss(pair);
            std::string token;
            std::vector<std::string> chr_cov;

            while (std::getline(ss, token, ':'))
            {
                chr_cov.push_back(token);
            }

            // Check if the pair is valid
            if (chr_cov.size() == 2)
            {
                // Get the chromosome and coverage
                std::string chr = chr_cov[0];
                double cov = std::stod(chr_cov[1]);

                // Add the pair to the map
                this->chr_cov[chr] = cov;

                // Print the pair
                std::cout << "Set mean coverage for " << chr << " to " << cov << std::endl;
            }
        }
    }
}

int InputData::getChrCov(std::string chr, double &cov)
{
    // Check if the chromosome is in the map
    if (this->chr_cov.find(chr) != this->chr_cov.end())
    {
        // Get the coverage
        cov = this->chr_cov[chr];

        // Return 0 if the chromosome is in the map
        return 0;
    }
    else
    {
        // Return -1 if the chromosome is not in the map
        return -1;
    }
}

std::string InputData::getPFBFilepath()
{
    return this->pfb_filepath;
}

void InputData::setPFBFilepath(std::string filepath)
{
    this->pfb_filepath = filepath;

    // Check if empty string
    if (filepath == "")
    {
        return;
        
    } else {
        // Check if the file exists
        FILE *fp = fopen(filepath.c_str(), "r");
        if (fp == NULL)
        {
            std::cerr << "PFB file does not exist: " << filepath << std::endl;
            exit(1);
        }
    }
}

void InputData::setThreadCount(int thread_count)
{
    this->thread_count = thread_count;
}

int InputData::getThreadCount()
{
    return this->thread_count;
}

std::string InputData::getHMMFilepath()
{
    return this->hmm_filepath;
}

void InputData::setHMMFilepath(std::string filepath)
{
    // Check if empty string
    if (filepath == "")
    {
        std::cout << "Using default HMM file: " << this->hmm_filepath << std::endl;
        return;
        
    } else {
        // Check if the file exists
        FILE *fp = fopen(filepath.c_str(), "r");
        if (fp == NULL)
        {
            std::cerr << "HMM file does not exist: " << filepath << std::endl;
            exit(1);
        } else {
            this->hmm_filepath = filepath;
            std::cout << "Using HMM file: " << this->hmm_filepath << std::endl;
        }
    }
}

void InputData::setDisableCIGAR(bool disable_cigar)
{
    this->disable_cigar = disable_cigar;
}

bool InputData::getDisableCIGAR()
{
    return this->disable_cigar;
}

void InputData::setCNVFilepath(std::string filepath)
{
    this->cnv_filepath = filepath;
}

std::string InputData::getCNVFilepath()
{
    return this->cnv_filepath;
}
