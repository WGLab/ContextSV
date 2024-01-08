#include "utils.h"

/// @cond
#include <stdio.h>
#include <string>
#include <iostream>
/// @endcond

#define BUFFER_SIZE 4096

// Print a progress bar
void printProgress(int progress, int total)
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

// Run bcftools to determine the chr notation of a VCF file
bool isChrNotation(std::string vcf_filepath)
{
    // Create the command to extract the chromosomes from the VCF
    std::string cmd = "bcftools query -f '%CHROM\n' " + vcf_filepath + " 2>/dev/null";

    // Open the pipe
    std::cout << "Checking if VCF is in chr notation..." << std::endl;
    std::cout << "Running command: " << cmd << std::endl;
    FILE* pipe = popen(cmd.c_str(), "r");
    if (pipe == NULL)
    {
        std::cerr << "Error: could not open pipe" << std::endl;
        return false;
    }

    // Read the first line
    char buffer[BUFFER_SIZE];
    if (!fgets(buffer, BUFFER_SIZE, pipe))
    {
        std::cerr << "Error reading from pipe" << std::endl;
        return false;
    }

    // Check if the first line contains "chr" using std::string::find
    bool is_chr_notation = false;
    std::string line(buffer);
    if (line.find("chr") != std::string::npos)
    {
        is_chr_notation = true;
    }
    
    // Close the pipe
    pclose(pipe);
    std::cout << "Done." << std::endl;

    return is_chr_notation;
}
