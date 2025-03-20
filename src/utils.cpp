#include "utils.h"

/// @cond
#include <sys/resource.h>  // getrusage
#include <iomanip>
#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <filesystem>
/// @endcond


// Define print mutex
std::mutex print_mtx;

// Print a progress bar
void printProgress(int progress, int total)
{
    float percent = (float)progress / (float)total * 100.0;
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
    FILE* pipe = popen(cmd.c_str(), "r");
    if (pipe == NULL)
    {
        std::cerr << "Error: could not open pipe" << std::endl;
        return false;
    }

    // Read the first line
    const int line_size = 256;
    char buffer[line_size];
    if (!fgets(buffer, line_size, pipe))
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

    return is_chr_notation;
}

// Thread-safe print message function
void printMessage(std::string message)
{
    std::lock_guard<std::mutex> lock(print_mtx);
    std::cout << message << std::endl;
}

// Thread-safe print error function
void printError(std::string message)
{
    std::lock_guard<std::mutex> lock(print_mtx);
    std::cerr << message << std::endl;
}

// Return the elapsed time given a start and end time (hours:minutes:seconds)
std::string getElapsedTime(std::chrono::high_resolution_clock::time_point start, std::chrono::high_resolution_clock::time_point end)
{
    std::chrono::duration<double> elapsed = end - start;
    int hours = elapsed.count() / 3600;
    int minutes = (elapsed.count() - (hours * 3600)) / 60;
    int seconds = elapsed.count() - (hours * 3600) - (minutes * 60);
    std::string elapsed_time = std::to_string(hours) + ":" + std::to_string(minutes) + ":" + std::to_string(seconds);
    return elapsed_time;
}

// Function to remove the 'chr' prefix from chromosome names
std::string removeChrPrefix(std::string chr)
{
    if (chr.find("chr") != std::string::npos) {
        return chr.substr(3);
    }
    return chr;
}

void printMemoryUsage(const std::string& functionName) {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);

    // Convert from KB to GB
    double mem_usage_gb = (double)usage.ru_maxrss / 1024.0 / 1024.0;
    std::cout << functionName << " memory usage: "
              << std::fixed << std::setprecision(2) << mem_usage_gb << " GB" << std::endl;
}

bool fileExists(const std::string &filepath)
{
    std::ifstream file(filepath);
    return file.is_open();
}

bool isFileEmpty(const std::string &filepath)
{
    return std::filesystem::file_size(filepath) == 0;
}

void closeJSON(const std::string &filepath)
{
    std::ofstream
        json_file(filepath, std::ios::app);

    json_file << "}\n";  // Close the last JSON object
    json_file << "]";  // Close the JSON array
    json_file.close();
}
