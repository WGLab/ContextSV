#include "sv_object.h"
#include "sv_object.h"
#include <algorithm>
#include <tuple>
#include <memory>
#include <cmath>
#include <stdexcept>
#include <iostream>

bool SVCall::operator<(const SVCall & other) const
{
    return std::tie(start, end, sv_type, alt_allele, data_type, genotype, hmm_likelihood) <
           std::tie(other.start, other.end, other.sv_type, other.alt_allele, other.data_type, other.genotype, other.hmm_likelihood);
}

void addSVCall(std::set<SVCall>& sv_calls, uint32_t start, uint32_t end, std::string sv_type, std::string alt_allele, std::string data_type, std::string genotype, double hmm_likelihood)
{
    // Throw an error if unknown SV type
    if (sv_type == "UNKNOWN") {
        throw std::runtime_error("ERROR: Cannot add unknown SV type");
    }
    
    if (start >= end) {
        throw std::runtime_error("ERROR: Invalid SV at position " + std::to_string(start) + "-" + std::to_string(end));
    }

    // If the SV call already exists (start and end position), then update all information if the
    // likelihood is higher
    // std::cout << "[TEST1] Adding SV call: " << start << "-" << end << " " << sv_type << " " << alt_allele << " " << data_type << " " << genotype << " " << hmm_likelihood << std::endl;
    std::vector<SVCall> updates;
    bool print_out = false;
    for (auto it = sv_calls.begin(); it != sv_calls.end();)
    {
        if (it->start == start && it->end == end)
        {
            if (hmm_likelihood > it->hmm_likelihood)
            {
                std::cout << "[DEBUG] Found higher likelihood for SV call: " << start << "-" << end << " " << sv_type << " " << alt_allele << " " << data_type << " " << genotype << " " << hmm_likelihood << std::endl;
                print_out = true;
                // Update the data type and support
                std::string new_data_type = it->data_type + "," + data_type;
                int new_support = it->support + 1;

                updates.push_back(SVCall{start, end, sv_type, alt_allele, new_data_type, genotype, hmm_likelihood, new_support});

                // Erase and re-insert the SV call
                // Erase the current iterator and safely insert the new SV call
                // sv_calls.erase(it);
                it = sv_calls.erase(it);  // Erase and get the next iterator
                // sv_calls.insert(SVCall{start, end, sv_type, alt_allele, new_data_type, genotype, hmm_likelihood, new_support});
            } else {
                // Return if no update is needed
                return;
            }
        } else {
            // Increment the iterator if the SV call does not match
            ++it;
        }
    }

    if (print_out)
    {
        std::cout << "[DEBUG] Adding updates" << std::endl;
    }
    
    // Insert the updates
    for (const auto& update : updates)
    {
        sv_calls.insert(update);
    }

    if (print_out)
    {
        std::cout << "[DEBUG] Added updates" << std::endl;
    }


    // Add the SV call if it does not exist
    // std::cout << "[TEST2] Adding SV call: " << start << "-" << end << " " << sv_type << " " << alt_allele << " " << data_type << " " << genotype << " " << hmm_likelihood << std::endl;
    // sv_calls.insert(SVCall{start, end, sv_type, alt_allele, data_type, genotype, hmm_likelihood, 1});
    // std::cout << "[TEST3] Added SV call: " << start << "-" << end << " " << sv_type << " " << alt_allele << " " << data_type << " " << genotype << " " << hmm_likelihood << std::endl;
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
