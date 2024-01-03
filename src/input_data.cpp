#include "input_data.h"

/// @cond
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <sstream>
#include <thread>

#include "utils.h"
/// @endcond

#define BUFFER_SIZE 4096
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
    this->disable_cigar = false;
    this->disable_snp_cnv = false;
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

        // If a region is set, load only the chromosome in the region (non-chr notation)
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

            std::cout << "Loading population allele frequency file for chromosome " << target_chr << std::endl;
        }

        // Load the PFB file and create the PFB map (chr -> VCF file)
        std::cout << "Loading population allele frequency file: " << filepath << std::endl;
        char buffer[BUFFER_SIZE];
        std::vector<std::thread> threads;  // Vector of threads for parallelization
        std::mutex pfb_mtx;  // Mutex for the PFB map
        std::mutex print_mtx;  // Mutex for printing messages
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
                // Get the chromosome
                std::string chr = chr_vcf[0];

                // If the region is set, skip chromosomes not in the region
                if (this->region_set && chr != target_chr)
                {
                    //std::cout << "Skipping chromosome " << chr << std::endl;
                    continue;
                } else {
                    std::cout << "Loading chromosome " << chr << std::endl;
                }

                // Get the VCF file
                std::string vcf = chr_vcf[1].substr(0, chr_vcf[1].find_first_of("\r\n"));  // Remove the newline character

                // Create a thread for processing the VCF file for each
                // chromosome
                threads.push_back(std::thread(&InputData::readChromosomeAFs, this, chr, vcf, std::ref(pfb_mtx), std::ref(print_mtx)));
            }
        }

        // Join the threads
        if (threads.size() == 0)
        {
            std::cerr << "Error: No chromosomes found in PFB file" << std::endl;
            exit(1);
        } else {
            std::cout << "Loading " << threads.size() << " chromosomes" << std::endl;
        }

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

void InputData::readChromosomeAFs(std::string chr, std::string filepath, std::mutex &pfb_mtx, std::mutex &print_mtx)
{
    // Check if the file exists
    this->printMessage("Checking that allele frequency file exists: " + filepath, print_mtx);
    FILE *fp = fopen(filepath.c_str(), "r");
    if (fp == NULL)
    {
        //std::cerr << "Error: Allele frequency file does not exist: " <<
        //filepath << std::endl;
        this->printError("Error: Allele frequency file does not exist: " + filepath, print_mtx);
        exit(1);
    }
    // Close the file
    fclose(fp);
    this->printMessage("Complete.", print_mtx);

    // Run bcftools index on the VCF file to create the index file if it does
    // not exist
    this->printMessage("Checking that allele frequency index file exists: " + filepath + ".csi", print_mtx);
    std::string index_cmd = "bcftools index -f " + filepath;
    system(index_cmd.c_str());
    this->printMessage("Complete.", print_mtx);

    // Check if the chromosome is in the reference genome and return it with the
    // reference notation
    this->printMessage("Checking that chromosome " + chr + " is in reference genome...", print_mtx);
    std::string chr_check = this->fasta_query.hasChromosome(chr);
    this->printMessage("Complete.", print_mtx);
    if (chr_check == "")
    {
        this->printError("Error: Chromosome " + chr + " not in reference genome", print_mtx);
        exit(1);
    }
    chr = chr_check;  // Update the chromosome with the reference notation

    // Load the allele frequency file and create the allele frequency map (position -> allele frequency)
    this->printMessage("Loading file: " + filepath, print_mtx);
    int af_count = 0;
    int af_min_hit = 0;  // Number of positions with allele frequency below the minimum
    int af_max_hit = 0;  // Number of positions with allele frequency above the maximum

    // Use bcftools to read the VCF file and create the allele frequency map
    // If a region is set, use the region
    std::string cmd;
    if (this->region == "")
    {
        cmd = "bcftools query -f '%POS\t%AF\n' -i 'INFO/variant_type=\"snv\"' " + filepath;
    } else {
        // Determine if the VCF file is in chr notation
        bool chr_notation = isChrNotation(filepath);

        // Fix the region's notation if necessary
        std::string target_region = this->region;
        if (chr_notation && this->region.find("chr") == std::string::npos)
        {
            // Add the chr notation
            target_region = "chr" + this->region;
        }
        else if (!chr_notation && this->region.find("chr") != std::string::npos)
        {
            // Remove the chr notation
            target_region = this->region.substr(3, this->region.size() - 3);
        }
        cmd = "bcftools query -f '%POS\t%AF\n' -i 'INFO/variant_type=\"snv\"' -r " + target_region + " " + filepath;
    }

    std::cout << "Command: " << cmd << std::endl;

    // Open a pipe to read the output of the command
    std::cout << "Opening pipe to read VCF file" << std::endl;
    FILE *pipe = popen(cmd.c_str(), "r");
    if (pipe == NULL)
    {
        this->printError("Error: Could not open pipe to read VCF file", print_mtx);
        exit(1);
    }
    std::cout << "Parsing AF values for chromosome " << chr << std::endl;

    // Create a PFB map for the chromosome
    std::map<int, double> chr_pfb_map;

    // Read the output of the command in a highly optimized manner
    char buffer[BUFFER_SIZE];
    std::string line;
    std::istringstream ss;
    std::string token;
    while (fgets(buffer, BUFFER_SIZE, pipe) != NULL)
    {
        // Remove the newline character
        line = buffer;
        line = line.substr(0, line.find_first_of("\r\n"));

        // Check if the line is valid
        if (line != "")
        {
            // Split the line by tab (position, allele frequency)
            ss.str(line);
            try
            {
                std::getline(ss, token, '\t');
                int pos = std::stoi(token);
                std::getline(ss, token, '\t');
                double af = std::stod(token);

                // Check if the allele frequency is within the valid range
                if (af >= MIN_PFB && af <= MAX_PFB)
                {
                    // Add the position and allele frequency to the map
                    chr_pfb_map[pos] = af;
                }
                else if (af < MIN_PFB)
                {
                    af_min_hit++;
                    chr_pfb_map[pos] = MIN_PFB;
                }
                else if (af > MAX_PFB)
                {
                    af_max_hit++;
                    chr_pfb_map[pos] = MAX_PFB;
                }
                af_count++;
            }
            catch (const std::invalid_argument &ia)
            {
                // Continue if the line is invalid. Usually this is due to the
                // allele frequency being empty '.' or 'NA'
            }

            // Clear the objects
            ss.clear();
            token.clear();
            line.clear();
        }
    }

    // Close the pipe
    pclose(pipe);
    std::cout << "Loaded " << af_count << " positions" << std::endl;

    // Check if the PFB map is empty
    if (chr_pfb_map.size() == 0)
    {
        this->printError("Error: No positions found for chromosome " + chr, print_mtx);
        exit(1);
    }

    // Add the chromosome PFB map to the PFB map
    this->addChromosomePopulationFrequency(chr, chr_pfb_map, pfb_mtx);

    // Print the number of positions found for the chromosome
    this->printMessage("Loaded " + std::to_string(this->pfb_map[chr].size()) + " positions for chromosome " + chr, print_mtx);

    // // Log the percentage of PFB values that were fixed
    // std::cout << "AF value count: " << af_count << std::endl;  // DEBUG
    // std::cout << "AF min. hit: " << af_min_hit << std::endl;  // DEBUG
    // std::cout << "AF max. hit: " << af_max_hit << std::endl;  // DEBUG
    // std::cout << "SNP AF values fixed: " << ((double) (af_min_hit + af_max_hit) / (double) af_count) * 100 << "%" << std::endl;
    // std::cout << "Min. fixed: " << ((double) af_min_hit / (double) af_count) * 100 << "%" << std::endl;
    // std::cout << "Max. fixed: " << ((double) af_max_hit / (double) af_count) * 100 << "%" << std::endl;
}

void InputData::addChromosomePopulationFrequency(std::string chr, std::map<int, double> pfb_map, std::mutex &mutex)
{
    std::lock_guard<std::mutex> lock(mutex);
    this->pfb_map[chr] = pfb_map;
}

void InputData::printMessage(std::string message, std::mutex &mutex)
{
    std::lock_guard<std::mutex> lock(mutex);
    std::cout << message << std::endl;
}

void InputData::printError(std::string message, std::mutex &mutex)
{
    std::lock_guard<std::mutex> lock(mutex);
    std::cerr << message << std::endl;
}
