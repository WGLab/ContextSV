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
	return start < other.start || (start == other.start && end < other.end);
    //return std::tie(start, end) < std::tie(other.start, other.end);
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

    // If the SV call already exists (start and end position), then update all information if the
    // likelihood is higher
    // std::cout << "[TEST1] Adding SV call: " << start << "-" << end << " " <<
    // sv_type << " " << alt_allele << " " << data_type << " " << genotype << "
    // " << hmm_likelihood << std::endl;
    sv_calls.insert(SVCall{start, end, sv_type, alt_allele, data_type, genotype, hmm_likelihood, 1});
    // SVCall new_sv_call{start, end, sv_type, alt_allele, data_type, genotype, hmm_likelihood, 1};
    
    // sv_calls.insert(new_sv_call);
    
    /*
    bool exists = false;
    bool print_out = false;
    for (auto it = sv_calls.begin(); it != sv_calls.end();)
    {
        if (it->start == start && it->end == end)
        {
            exists = true;
            if (hmm_likelihood > it->hmm_likelihood)
            {
                //std::cout << "[DEBUG] Found higher likelihood for SV call: " << start << "-" << end << " " << sv_type << " " << alt_allele << " " << data_type << " " << genotype << " " << hmm_likelihood << std::endl;
                print_out = true;
                // Update the data type and support
                // std::string new_data_type = it->data_type + "," + data_type;
                // int new_support = it->support + 1;
                new_sv_call.data_type = it->data_type + "," + data_type;
                new_sv_call.support = it->support + 1;
                //higher_lh = true;

                // updates.push_back(SVCall{start, end, sv_type, alt_allele, new_data_type, genotype, hmm_likelihood, new_support});

                // Erase and re-insert the SV call
                // Erase the current iterator and safely insert the new SV calls
                std::cout << "Erasing iterator." << std::endl;
                sv_calls.erase(it);
                std::cout << "Iterator erased." << std::endl;
                break;
                //it = sv_calls.erase(it);  // Erase and get the next iterator
                // sv_calls.insert(SVCall{start, end, sv_type, alt_allele, new_data_type, genotype, hmm_likelihood, new_support});
            } else {
                // End if the SV exists but is lower lh
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

    // Update the SV call if it does not exist, or if the likelihood is higher
    // than the existing call
    if (print_out)
    {
        std::cout << "[DEBUG] Inserting call" << std::endl;
    }
    sv_calls.insert(new_sv_call);
    if (print_out)
    {
        std::cout << "[DEBUG] Call inserted" << std::endl;
    }
    // Insert the updates
    // for (const auto& update : updates)
    // {
    //     sv_calls.insert(update);
    // }

    // if (print_out)
    // {
    //     std::cout << "[DEBUG] Added updates" << std::endl;
    // }


    // Add the SV call if it does not exist
    // std::cout << "[TEST2] Adding SV call: " << start << "-" << end << " " << sv_type << " " << alt_allele << " " << data_type << " " << genotype << " " << hmm_likelihood << std::endl;
    // sv_calls.insert(SVCall{start, end, sv_type, alt_allele, data_type, genotype, hmm_likelihood, 1});
    // std::cout << "[TEST3] Added SV call: " << start << "-" << end << " " << sv_type << " " << alt_allele << " " << data_type << " " << genotype << " " << hmm_likelihood << std::endl;
    */
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
    // int initial_size = sv_calls.size();
    std::vector<SVCall> merged_sv_calls;
    auto it = sv_calls.begin();
    SVCall current_merge = *it++;

    for (; it != sv_calls.end(); ++it) {
        const SVCall& next = *it;

        // Check if the SV calls overlap by at least 50%
        uint32_t overlap_start = std::max(current_merge.start, next.start);
        uint32_t overlap_end = std::min(current_merge.end, next.end);
        uint32_t overlap_length = (overlap_start < overlap_end) ? overlap_end - overlap_start : 0;

        uint32_t current_length = current_merge.end - current_merge.start;
        uint32_t next_length = next.end - next.start;

        // Merge the SV calls if the overlap is > 0
        //double overlap_pct_current = static_cast<double>(overlap_length) / current_length;
        //double overlap_pct_next = static_cast<double>(overlap_length) / next_length;

        //if (overlap_pct_current >= 0.5 || overlap_pct_next >= 0.5) {
        if (overlap_length > 0) {
            // Merge the SV calls based on the likelihood
            if (next.hmm_likelihood != 0.0) {
                // Update the likelihood if the next SV call has a likelihood
                // and it is higher than the current merged SV call
                if (next.hmm_likelihood > current_merge.hmm_likelihood) {
                    current_merge = next;
                }
            } else {
                // If both have no likelihood (CIGAR only), then merge the SV calls
                // based on largest SV length
                if (next.hmm_likelihood == current_merge.hmm_likelihood) {
                    if (next_length > current_length) {
                        current_merge = next;
                    }
                }
                // if (next_length > current_length) {
                //     current_merge = next;
                // }
            }
        } else {
            // No overlap: Save the previous SV and continue
            merged_sv_calls.push_back(current_merge);
            current_merge = next;
        }
    }

    // Add the last merged SV call
    merged_sv_calls.push_back(current_merge);

    // Update the SV calls
    sv_calls.clear();
    for (const auto& sv_call : merged_sv_calls) {
        sv_calls.insert(sv_call);
    }
    // int updated_size = sv_calls.size();
    // std::cout << "Merged " << initial_size << " SV calls into " << updated_size << " SV calls" << std::endl;
}
