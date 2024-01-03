// Utility functions

#ifndef UTILS_H
#define UTILS_H

/// @cond
#include <string>
/// @endcond


// Print the progress of a task
void printProgress(int progress, int total);

// Run bcftools to determine the chr notation of a VCF file
bool isChrNotation(std::string vcf_filepath);

#endif // UTILS_H
