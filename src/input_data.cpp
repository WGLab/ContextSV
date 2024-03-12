#include "input_data.h"

/// @cond
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <sstream>
#include <thread>

#include "utils.h"
/// @endcond

#define MIN_PFB 0.01  // Minimum SNP population allele frequency
#define MAX_PFB 0.99  // Maximum SNP population allele frequency

// Constructor
InputData::InputData()
{
    this->short_read_bam = "";
    this->long_read_bam = "";
    this->ref_filepath = "";
    this->snp_vcf_filepath = "";
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
    this->disable_cigar = false;
    this->disable_snp_cnv = false;
    this->verbose = false;
    this->save_cnv_data = false;
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

const FASTAQuery &InputData::getRefGenome() const
{
    return this->fasta_query;
}

std::string InputData::queryRefGenome(std::string chr, int64_t pos_start, int64_t pos_end)
{
    return this->fasta_query.query(chr, pos_start, pos_end);
}

std::vector<std::string> InputData::getRefGenomeChromosomes()
{
    return this->fasta_query.getChromosomes();
}

int64_t InputData::getRefGenomeChromosomeLength(std::string chr)
{
    return this->fasta_query.getChromosomeLength(chr);
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
        this->region_set = true;

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

void InputData::setMeanChromosomeCoverage(std::string chr_cov)
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

double InputData::getMeanChromosomeCoverage(std::string chr)
{
    // Using find to check if the key exists
    auto it = chr_cov.find(chr);

    // If key is not found, throw an error
    if (it == chr_cov.end()) {
        throw std::out_of_range("Key not found in the map.");
    }

    // Key exists, return the corresponding double value
    return it->second;
}

void InputData::setAlleleFreqFilepaths(std::string filepath)
{
    // this->pfb_filepath = filepath;

    // Check if empty string
    if (filepath == "")
    {
        return;
        
    } else {
        // Check if the file exists
        FILE *fp = fopen(filepath.c_str(), "r");
        if (fp == NULL)
        {
            std::cerr << "Population allele frequency file does not exist: " << filepath << std::endl;
            exit(1);
        }

        // If a region is set, load only the chromosome in the region
        std::string target_chr;
        if (this->region_set)
        {
            target_chr = this->region_chr;

            // Check if the region is in chr notation
            if (target_chr.find("chr") != std::string::npos)
            {
                // Remove the chr notation
                target_chr = target_chr.substr(3, target_chr.size() - 3);
            }
            //std::cout << "Loading population allele frequency file for chromosome " << target_chr << std::endl;
        }

        // Load the file and create a map of chromosome -> VCF file
        if (this->verbose)
        {
            std::cout << "Loading population allele frequency files from: " << filepath << std::endl;
        }

        // Read the file line by line
        const int line_size = 256;  // Sufficient buffer size for each line
        char buffer[line_size];
        while (fgets(buffer, line_size, fp) != NULL)
        {
            // Check if the line is valid
            if (buffer[0] == '#')
            {
                continue;
            }

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
                // Get the chromosome
                std::string chr = chr_vcf[0];

                // Get the VCF file
                std::string vcf = chr_vcf[1].substr(0, chr_vcf[1].find_first_of("\r\n"));  // Remove the newline character

                // Check if the file exists
                FILE *vcf_fp = fopen(vcf.c_str(), "r");
                if (vcf_fp == NULL)
                {
                    std::cerr << "Error: Allele frequency file does not exist: " << vcf << std::endl;
                    exit(1);
                }
                // Close the file
                fclose(vcf_fp);

                // Add the chromosome and VCF file to the map
                this->pfb_filepaths[chr] = vcf;
            }
        }
    }
}

std::string InputData::getAlleleFreqFilepath(std::string chr)
{
    // Remove the chr notation
    if (chr.find("chr") != std::string::npos)
    {
        chr = chr.substr(3, chr.size() - 3);
    }
    return this->pfb_filepaths[chr];
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

void InputData::setDisableSNPCNV(bool disable_snp_cnv)
{
    this->disable_snp_cnv = disable_snp_cnv;
}

bool InputData::getDisableSNPCNV()
{
    return this->disable_snp_cnv;
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

void InputData::setVerbose(bool verbose)
{
    this->verbose = verbose;
}

bool InputData::getVerbose()
{
    return this->verbose;
}

void InputData::saveCNVData(bool save_cnv_data)
{
    this->save_cnv_data = save_cnv_data;
}

bool InputData::getSaveCNVData()
{
    return this->save_cnv_data;
}
