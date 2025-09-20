// Utility functions

#ifndef UTILS_H
#define UTILS_H

#include <htslib/sam.h>
#include <htslib/synced_bcf_reader.h>

/// @cond
#include <string>
#include <mutex>
#include <chrono>
/// @endcond


// Print a message to stdout in a thread-safe manner
void printMessage(std::string message);

// Print an error message to stderr in a thread-safe manner
void printError(std::string message);

std::string getElapsedTime(std::chrono::high_resolution_clock::time_point start, std::chrono::high_resolution_clock::time_point end);

void printMemoryUsage(const std::string &functionName);

bool fileExists(const std::string &filepath);

bool isFileEmpty(const std::string &filepath);

void closeJSON(const std::string & filepath);

#endif // UTILS_H
