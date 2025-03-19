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


// Guard to close the BAM file
struct BamFileGuard {
    samFile* fp_in;
    hts_idx_t* idx;
    bam_hdr_t* bamHdr;

    BamFileGuard(samFile* fp_in, hts_idx_t* idx, bam_hdr_t* bamHdr)
        : fp_in(fp_in), idx(idx), bamHdr(bamHdr) {}

    ~BamFileGuard() {
        if (idx) {
            hts_idx_destroy(idx);
            idx = nullptr;
        }
        if (bamHdr) {
            bam_hdr_destroy(bamHdr);
            bamHdr = nullptr;
        }
        if (fp_in) {
            sam_close(fp_in);
            fp_in = nullptr;
        }
    }

    BamFileGuard(const BamFileGuard&) = delete;  // Non-copyable
    BamFileGuard& operator=(const BamFileGuard&) = delete;  // Non-assignable
};

// Print the progress of a task
void printProgress(int progress, int total);

// Run bcftools to determine the chr notation of a VCF file
bool isChrNotation(std::string vcf_filepath);

// Print a message to stdout in a thread-safe manner
void printMessage(std::string message);

// Print an error message to stderr in a thread-safe manner
void printError(std::string message);

std::string getElapsedTime(std::chrono::high_resolution_clock::time_point start, std::chrono::high_resolution_clock::time_point end);

std::string removeChrPrefix(std::string chr);

void printMemoryUsage(const std::string &functionName);

bool fileExists(const std::string &filepath);

void openJSON(const std::string & filepath);

void closeJSON(const std::string & filepath);

#endif // UTILS_H
