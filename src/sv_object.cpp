#include "sv_object.h"

#include <algorithm>
#include <tuple>
#include <memory>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <numeric>

#include "utils.h"

bool SVCall::operator<(const SVCall & other) const
{
	return start < other.start || (start == other.start && end < other.end);
}

void addSVCall(std::vector<SVCall>& sv_calls, uint32_t start, uint32_t end, SVType sv_type, const std::string& alt_allele, std::string data_type, std::string genotype, double hmm_likelihood, int read_depth, uint8_t qual)
{
    // Ignore unknown SV types
    // if (sv_type == "UNKNOWN" || sv_type == "NEUTRAL") {
    //     return;
    // }
    if (sv_type == SVType::UNKNOWN || sv_type == SVType::NEUTRAL) {
        return;
    }
    
    if (start > end) {
        printError("ERROR: Invalid SV at position " + std::to_string(start) + "-" + std::to_string(end));
        return;
    }

    // Insert the SV call in sorted order
    // SVCall sv_call{start, end, sv_type, alt_allele, data_type, genotype,
    // hmm_likelihood, read_depth, 1};
    SVCall sv_call{start, end, sv_type, alt_allele, data_type, genotype, hmm_likelihood, read_depth, 1, 1, qual};
    auto it = std::lower_bound(sv_calls.begin(), sv_calls.end(), sv_call);

    // Determine if the SV call already exists
    if (it != sv_calls.end() && it->start == start && it->end == end)
    {
        it->support += 1;  // Update the read support

        // Update SV type if likelihood is higher
        if (hmm_likelihood != 0.0 && hmm_likelihood > it->hmm_likelihood)
        {
            // Update the SV call
            it->sv_type = sv_type;
            it->data_type = data_type;
            it->genotype = genotype;
            it->hmm_likelihood = hmm_likelihood;
            it->qual = qual;
        }
    } else {
        sv_calls.insert(it, sv_call);  // Insert the new SV call
    }
}

uint32_t getSVCount(const std::vector<SVCall>& sv_calls)
{
    return (uint32_t) sv_calls.size();
}

void concatenateSVCalls(std::vector<SVCall> &target, const std::vector<SVCall>& source)
{
    target.insert(target.end(), source.begin(), source.end());
}

void mergeSVs(std::vector<SVCall>& sv_calls)
{
    if (sv_calls.size() < 2) {
        return;
    }
    int initial_size = sv_calls.size();

    std::vector<bool> merged(sv_calls.size(), false);
    std::vector<SVCall> merged_sv_calls;

    // Sort SVs by start position to improve efficiency
    std::sort(sv_calls.begin(), sv_calls.end(), [](const SVCall& a, const SVCall& b) {
        return a.start < b.start;
    });

    for (size_t i = 0; i < sv_calls.size(); i++) {
        if (merged[i]) continue;

        size_t best_index = i;
        int total_cluster_size = sv_calls[i].cluster_size;  // Track total cluster size

        for (size_t j = i + 1; j < sv_calls.size(); j++) {
            if (merged[j]) continue;

            // Compute overlap
            uint32_t overlap_start = std::max(sv_calls[i].start, sv_calls[j].start);
            uint32_t overlap_end = std::min(sv_calls[i].end, sv_calls[j].end);
            uint32_t overlap_length = (overlap_end > overlap_start) ? (overlap_end - overlap_start) : 0;

            // Compute union length correctly
            uint32_t union_start = std::min(sv_calls[i].start, sv_calls[j].start);
            uint32_t union_end = std::max(sv_calls[i].end, sv_calls[j].end);
            uint32_t union_length = union_end - union_start;  // No +1 to prevent off-by-one errors

            double overlap_fraction = (union_length > 0) ? (static_cast<double>(overlap_length) / union_length) : 0.0;

            // Throw error if fraction > 1
            if (overlap_fraction > 1.0) {
                throw std::runtime_error("Error: Overlap fraction = " + std::to_string(overlap_fraction) + " > 1.0");
            }

            // if (overlap_fraction > 0.5) {
            if (overlap_fraction > 0.5) {  // Changed from 0.5
                total_cluster_size += sv_calls[j].cluster_size;
                if (sv_calls[j].support > sv_calls[best_index].support) {
                    best_index = j;
                }
                merged[j] = true;  // Mark SV as merged
            }
        }

        sv_calls[best_index].cluster_size = total_cluster_size; // Update best SV with total size
        merged_sv_calls.push_back(sv_calls[best_index]); // Keep the strongest SV
    }

    // Filter out merged SVs with low support or cluster size
    // merged_sv_calls.erase(std::remove_if(merged_sv_calls.begin(), merged_sv_calls.end(), [initial_size](const SVCall& sv_call) {
    //     return sv_call.support < 2 && sv_call.cluster_size < 10;  // Adjust thresholds as needed
    // }), merged_sv_calls.end());
    // merged_sv_calls.erase(std::remove_if(merged_sv_calls.begin(), merged_sv_calls.end(), [initial_size](const SVCall& sv_call) {
    //     return sv_call.support < 2 && sv_call.cluster_size < 3;  // Adjust thresholds as needed
    // }), merged_sv_calls.end());

    sv_calls = std::move(merged_sv_calls); // Replace with filtered list

    // Print SVs that have length 2039
    for (const auto& sv_call : sv_calls) {
        if (sv_call.end - sv_call.start == 2039) {
            printMessage("Found merged SV with length 2039 at " + std::to_string(sv_call.start) + "-" + std::to_string(sv_call.end) + " (SUP=" + std::to_string(sv_call.support) + ")");
        }
    }

    int updated_size = sv_calls.size();
    printMessage("Merged " + std::to_string(initial_size) + " SV calls into " + std::to_string(updated_size) + " SV calls");
}

void filterSVsWithLowSupport(std::vector<SVCall>& sv_calls, int min_support)
{
    // Filter SV calls with low read support or low cluster size
    sv_calls.erase(std::remove_if(sv_calls.begin(), sv_calls.end(), [min_support](const SVCall& sv_call) {
        return sv_call.support < min_support && sv_call.cluster_size < min_support;
    }), sv_calls.end());
}

void filterSVsWithLowSupport(std::vector<SVCall> &sv_calls, int min_support, const std::string &data_type)
{
    // Filter SV calls with low read depth only for the specified data type, keeping the rest
    sv_calls.erase(std::remove_if(sv_calls.begin(), sv_calls.end(), [min_support, data_type](const SVCall& sv_call) {
        return sv_call.support < min_support && sv_call.data_type == data_type;
    }), sv_calls.end());
}
