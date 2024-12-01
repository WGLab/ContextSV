#include "sv_object.h"

#include <algorithm>
#include <tuple>
#include <memory>
#include <cmath>
#include <stdexcept>
#include <iostream>

#include "utils.h"

bool SVCall::operator<(const SVCall & other) const
{
	return start < other.start || (start == other.start && end < other.end);
}

void addSVCall(std::set<SVCall>& sv_calls, uint32_t start, uint32_t end, std::string sv_type, std::string alt_allele, std::string data_type, std::string genotype, double hmm_likelihood)
{
    // Ignore unknown SV types
    if (sv_type == "UNKNOWN") {
        return;
    }
    
    if (start >= end) {
        throw std::runtime_error("ERROR: Invalid SV at position " + std::to_string(start) + "-" + std::to_string(end));
    }

    // printMessage("Adding SV call: " + std::to_string(start) + "-" + std::to_string(end) + " with length " + std::to_string(end - start) + " and type " + sv_type);
    sv_calls.insert(SVCall{start, end, sv_type, alt_allele, data_type, genotype, hmm_likelihood, 1});
}

std::vector<std::set<SVCall>> splitSVsIntoChunks(std::set<SVCall>& sv_calls, int chunk_count)
{
    // Split the SV calls into chunks
    std::vector<std::set<SVCall>> sv_chunks;
    int sv_count = (int) sv_calls.size();
    int chunk_size = std::ceil((double) sv_count / (double) chunk_count);
    int current_chunk = 0;
    std::set<SVCall> current_sv_chunk;
    for (const auto& sv_call : sv_calls)
    {
        current_sv_chunk.insert(sv_call);

        // If the current chunk size is reached, then add the chunk to the
        // vector and reset the current chunk
        if ((int) current_sv_chunk.size() == chunk_size)
        {
            // sv_chunks.insert(current_sv_chunk);
            sv_chunks.push_back(current_sv_chunk);
            current_sv_chunk.clear();
            current_chunk++;
        }
    }

    // Add the last chunk if it is not empty
    if (!current_sv_chunk.empty())
    {
        sv_chunks.push_back(current_sv_chunk);
        // sv_chunks.insert(current_sv_chunk);
    }

    return sv_chunks;
}

uint32_t getSVCount(const std::set<SVCall>& sv_calls)
{
    return (uint32_t) sv_calls.size();
}

void concatenateSVCalls(std::set<SVCall> &target, const std::set<SVCall> &source)
{
    // Efficiently concatenate two sets of SV calls
    target.insert(source.begin(), source.end());
}

void mergeSVs(std::set<SVCall>& sv_calls) {
    if (sv_calls.size() < 2) {
        return;
    }

    // Merge SV calls if they overlap by at least 50%
    int initial_size = sv_calls.size();
    std::vector<SVCall> merged_sv_calls;
    auto it = sv_calls.begin();
    SVCall current_merge = *it++;
    for (; it != sv_calls.end(); ++it) {
        const SVCall& next = *it;

        // Find overlap
        if (next.start <= current_merge.end) {
            // Merge the SV calls if it is a subset
            if (next.end <= current_merge.end) {
                continue;
            }

            // Merge the SV calls based on HMM log likelihood (keep the higher
            // likelihood), 0.0 indicates no likelihood
            if (next.hmm_likelihood != 0.0 && next.hmm_likelihood > current_merge.hmm_likelihood) {
                current_merge = next;  // Continue with the next call
            }

        } else {
            // No overlap: Save the previous SV and continue
            merged_sv_calls.push_back(current_merge);
            current_merge = next;
        }
    }

    // Add the last merged SV call
    // printMessage("Saving SV call: " + std::to_string(current_merge.start) + "-" + std::to_string(current_merge.end) + " with likelihood " + std::to_string(current_merge.hmm_likelihood));
    merged_sv_calls.push_back(current_merge);

    // Replace contents of the SV calls
    sv_calls = std::set<SVCall>(merged_sv_calls.begin(), merged_sv_calls.end());

    // // Update the SV calls
    // sv_calls.clear();
    // for (const auto& sv_call : merged_sv_calls) {
    //     sv_calls.insert(sv_call);
    // }
    int updated_size = sv_calls.size();
    std::cout << "Merged " << initial_size << " SV calls into " << updated_size << " SV calls" << std::endl;
}
