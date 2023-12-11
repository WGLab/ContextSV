#include "input_data.h"

/// @cond
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <sstream>
#include <thread>
/// @endcond

#define BUFFER_SIZE 1024
#define MIN_PFB 0.01  // Minimum SNP population allele frequency
#define MAX_PFB 0.99  // Maximum SNP population allele frequency

// Constructor
InputData::InputData()
{
    this->short_read_bam = "";
    this->long_read_bam = "";
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
    this->whole_genome = true;
}

std::string InputData::getShortReadBam()
{
    return this->short_read_bam;
}

void InputData::setShortReadBam(std::string filepath)
{
    this->short_read_bam = filepath;

    // Check if empty string
    if (filepath == "")
    {
        return;
        
    } else {
        // Check if the file exists
        FILE *fp = fopen(filepath.c_str(), "r");
        if (fp == NULL)
        {
            std::cerr << "Short read BAM file does not exist: " << filepath << std::endl;
            exit(1);
        }
    }
}

std::string InputData::getLongReadBam()
{
    return this->long_read_bam;
}

void InputData::setLongReadBam(std::string filepath)
{
    this->long_read_bam = filepath;

    // Check if empty string
    if (filepath == "")
    {
        return;
        
    } else {
        // Check if the file exists
        FILE *fp = fopen(filepath.c_str(), "r");
        if (fp == NULL)
        {
            std::cerr << "Long read BAM file does not exist: " << filepath << std::endl;
            exit(1);
        }
    }
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
    
    // Check if empty string. This means all chromosomes will be used
    if (region == "")
    {
        std::cout << "No region specified. Using all chromosomes" << std::endl;
        this->whole_genome = true;
        return;
    }

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
        std::cerr << "Error: Could not parse region" << std::endl;
        exit(1);
    }

    // Check if only a chromosome was parsed
    if (this->region_start == 0 && this->region_end == 0)
    {
        // Use the entire chromosome as the region
        this->region_set = false;

        // Set the whole genome flag to false
        this->whole_genome = false;
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

            // Set the whole genome flag to false
            this->whole_genome = false;
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

void InputData::setAlleleFreqFilepaths(std::string filepath)
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

        // Load the PFB file and create the PFB map (chr -> VCF file)
        std::cout << "Loading population allele frequency file: " << filepath << std::endl;
        char buffer[BUFFER_SIZE];
        std::vector<std::thread> threads;  // Vector of threads for parallelization
        while (fgets(buffer, BUFFER_SIZE, fp) != NULL)
        {
            // Split by equals sign (chromosome, VCF file)
            std::istringstream ss(buffer);
            std::string token;
            std::vector<std::string> chr_vcf;

            while (std::getline(ss, token, '='))
            {
                chr_vcf.push_back(token);
            }

            // Check if the line is valid
            if (chr_vcf.size() == 2)
            {
                // Get the chromosome and VCF file
                std::string chr = chr_vcf[0];
                std::string vcf = chr_vcf[1];

                // Remove the newline and null characters
                vcf[strcspn(vcf.c_str(), "\r")] = 0;
                vcf[strcspn(vcf.c_str(), "\n")] = 0;

                // Read the VCF file and create the allele frequency map
                // (position -> allele frequency)

                // Read the VCF file and create the allele frequency map

                // Create a thread for processing the VCF file for each chromosome
                threads.push_back(std::thread(&InputData::readChromosomeAFs, this, chr, vcf));
            }
        }

        // Join the threads
        for (auto &thread : threads)
        {
            thread.join();
        }

        // Close the file
        fclose(fp);
        std::cout << "Loaded " << this->pfb_map.size() << " chromosomes" << std::endl;
    }
}

std::string InputData::getAlleleFreqFilepaths()
{
    return this->pfb_filepath;
}

PFBMap InputData::getPFBMap()
{
    return this->pfb_map;
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

void InputData::setWholeGenome(bool whole_genome)
{
    this->whole_genome = whole_genome;
}

bool InputData::getWholeGenome()
{
    return this->whole_genome;
}

void InputData::readChromosomeAFs(std::string chr, std::string filepath)
{
    // Check if the file exists
    FILE *fp = fopen(filepath.c_str(), "r");
    if (fp == NULL)
    {
        std::cerr << "Error: Allele frequency file does not exist: " << filepath << std::endl;
        exit(1);
    }

    // Check if the chromosome is in the reference genome and return it with the
    // reference notation    
    std::string chr_check = this->fasta_query.hasChromosome(chr);
    if (chr_check == "")
    {
        std::cerr << "Error: Chromosome " << chr << " not in reference genome" << std::endl;
        exit(1);
    }
    chr = chr_check;  // Update the chromosome with the reference notation

    // Load the allele frequency file and create the allele frequency map (position -> allele frequency)
    std::cout << "Loading allele frequency file: " << filepath << std::endl;
    char buffer[BUFFER_SIZE];
    int af_count = 0;
    int af_min_hit = 0;  // Number of positions with allele frequency below the minimum
    int af_max_hit = 0;  // Number of positions with allele frequency above the maximum
    while (fgets(buffer, BUFFER_SIZE, fp) != NULL)
    {
        // Remove the newline character
        buffer[strcspn(buffer, "\n")] = 0;

        // Split the line by tab (position, allele frequency)
        std::istringstream ss(buffer);
        std::string token;
        std::vector<std::string> pos_af;

        while (std::getline(ss, token, '\t'))
        {
            pos_af.push_back(token);
        }

        // Check if the line is valid
        if (pos_af.size() == 2)
        {
            // Get the position and allele frequency
            int pos = atoi(pos_af[0].c_str());
            double af = atof(pos_af[1].c_str());

            // Check if the allele frequency is within the valid range
            if (af >= MIN_PFB && af <= MAX_PFB)
            {
                // Add the position and allele frequency to the map
                this->pfb_map[chr][pos] = af;
                af_count++;
            }
            else if (af < MIN_PFB)
            {
                af_min_hit++;
            }
            else if (af > MAX_PFB)
            {
                af_max_hit++;
            }
        }
    }

    // Close the file
    fclose(fp);

    // Log the percentage of PFB values that were fixed
    std::cout << "SNP AF values fixed: " << (double) (af_min_hit + af_max_hit) / (double) af_count * 100 << "%" << std::endl;
    std::cout << "Min. fixed: " << (double) af_min_hit / (double) af_count * 100 << "%" << std::endl;
    std::cout << "Max. fixed: " << (double) af_max_hit / (double) af_count * 100 << "%" << std::endl;
}
