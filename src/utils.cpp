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

#ifndef VERSION
#define VERSION "vUNKNOWN"
#endif


std::mutex print_mtx;


static std::string run_cmd(const char* cmd) {
    std::array<char, 256> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) return {};
    while (fgets(buffer.data(), buffer.size(), pipe.get())) result += buffer.data();
    if (!result.empty() && result.back() == '\n') result.pop_back();
    return result;
}

std::string currentVersion() {
    auto gitv = run_cmd("git describe --tags --always 2>/dev/null");
    if (!gitv.empty()) return gitv;
    return VERSION;
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


