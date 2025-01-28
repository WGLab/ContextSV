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
    this->chr = "";
    this->start_end = std::make_pair(0, 0);
    this->region_set = false;
    this->output_dir = "";
    this->sample_size = 100;
    this->min_cnv_length = 1000;
    this->min_reads = 5;
    this->thread_count = 1;
    this->hmm_filepath = "data/wgs.hmm";
    this->verbose = false;
    this->save_cnv_data = false;
    this->single_chr = false;
}

std::string InputData::getShortReadBam() const
{
    return this->short_read_bam;
}

void InputData::setShortReadBam(std::string filepath)
{
    this->short_read_bam = filepath;

    // Check if empty string
    if (filepath.empty())
    {
        return;
        
    } else {
        // Check if the file exists
        FILE *fp = fopen(filepath.c_str(), "r");
        if (fp == NULL)
        {
            throw std::runtime_error("Short read BAM file does not exist: " + filepath);
        } else {
            fclose(fp);
        }
    }
}

std::string InputData::getLongReadBam() const
{
    return this->long_read_bam;
}

void InputData::setLongReadBam(std::string filepath)
{
    this->long_read_bam = filepath;

    // Check if empty string
    if (filepath.empty())
    {
        return;
        
    } else {
        // Check if the file exists
        FILE *fp = fopen(filepath.c_str(), "r");
        if (fp == NULL)
        {
            throw std::runtime_error("Long read BAM file does not exist: " + filepath);
        } else {
            fclose(fp);
        }
    }
}

void InputData::setRefGenome(std::string filepath)
{
    this->ref_filepath = filepath;
}

std::string InputData::getRefGenome() const
{
    return this->ref_filepath;
}

std::string InputData::getOutputDir() const
{
    return this->output_dir;
}

void InputData::setOutputDir(std::string dirpath)
{
    this->output_dir = dirpath;
    std::string cmd = "mkdir -p " + output_dir;
    try
    {
        std::system(cmd.c_str());
    } catch (const std::exception& e)
    {
        std::cerr << "Error creating output directory: " << e.what() << std::endl;
        exit(1);
    }
}

int InputData::getSampleSize() const
{
    return this->sample_size;
}

void InputData::setSampleSize(int sample_size)
{
    this->sample_size = sample_size;
}

std::string InputData::getSNPFilepath() const
{
    return this->snp_vcf_filepath;
}

void InputData::setSNPFilepath(std::string filepath)
{
    this->snp_vcf_filepath = filepath;
}

std::string InputData::getEthnicity() const
{
    return this->ethnicity;
}

void InputData::setEthnicity(std::string ethnicity)
{
    this->ethnicity = ethnicity;
}

uint32_t InputData::getMinCNVLength() const
{
    return this->min_cnv_length;
}

void InputData::setMinCNVLength(int min_cnv_length)
{
    this->min_cnv_length = (uint32_t) min_cnv_length;
}

void InputData::setMinReadSupport(int min_reads)
{
    // Ensure that the minimum read support is an integer and greater than 0
    if (min_reads < 1)
    {
        throw std::runtime_error("Minimum read support must be an integer greater than 0");
    }
    this->min_reads = min_reads;
}

int InputData::getMinReadSupport() const
{
    return this->min_reads;
}

void InputData::setChromosome(std::string chr)
{
    this->chr = chr;
    this->single_chr = true;
}

std::string InputData::getChromosome() const
{
    return this->chr;
}

bool InputData::isSingleChr() const
{
    return this->single_chr;
}

void InputData::setRegion(std::string region)
{
    // Check if the region is valid
    if (region != "")
    {
        // Split the region by colon
        std::istringstream ss(region);
        std::string token;
        std::vector<std::string> region_tokens;

        while (std::getline(ss, token, '-'))
        {
            region_tokens.push_back(token);
        }

        // Check if the region is valid
        if (region_tokens.size() == 2)
        {
            // Get the start and end positions
            int32_t start = std::stoi(region_tokens[0]);
            int32_t end = std::stoi(region_tokens[1]);

            // Set the region
            this->start_end = std::make_pair(start, end);
            this->region_set = true;

            std::cout << "Region set to " << this->chr << ":" << start << "-" << end << std::endl;
        }
    }
}

std::pair<int32_t, int32_t> InputData::getRegion() const
{
    return this->start_end;
}

bool InputData::isRegionSet() const
{
    return this->region_set;
}

void InputData::setAlleleFreqFilepaths(std::string filepath)
{
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
        if (this->chr != "")
        {
            target_chr = this->chr;

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

std::string InputData::getAlleleFreqFilepath(std::string chr) const
{
    // Remove the chr notation
    if (chr.find("chr") != std::string::npos)
    {
        chr = chr.substr(3, chr.size() - 3);
    }

    try
    {
        return this->pfb_filepaths.at(chr);
    }
    catch (const std::out_of_range& e)
    {
        return "";
    }
}

void InputData::setThreadCount(int thread_count)
{
    this->thread_count = thread_count;
}

int InputData::getThreadCount() const
{
    return this->thread_count;
}

std::string InputData::getHMMFilepath() const
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

bool InputData::getSaveCNVData() const
{
    return this->save_cnv_data;
}
