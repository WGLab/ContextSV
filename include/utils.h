// Utility functions

#ifndef UTILS_H
#define UTILS_H

/// @cond
#include <string>
#include <mutex>
/// @endcond


// Print the progress of a task
void printProgress(int progress, int total);

// Run bcftools to determine the chr notation of a VCF file
bool isChrNotation(std::string vcf_filepath);

// Print a message to stdout in a thread-safe manner
void printMessage(std::string message, std::mutex & mutex);

// Print an error message to stderr in a thread-safe manner
void printError(std::string message, std::mutex & mutex);

#endif // UTILS_H
