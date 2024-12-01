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

void addSVCall(std::vector<SVCall>& sv_calls, uint32_t start, uint32_t end, std::string sv_type, std::string alt_allele, std::string data_type, std::string genotype, double hmm_likelihood)
{
    // Catch underflow errors
    if (start > 4000000000 || end > 4000000000) {
        throw std::runtime_error("ERROR: Integer underflow for SV call at position " + std::to_string(start) + "-" + std::to_string(end));
    }

    // Ignore unknown SV types
    if (sv_type == "UNKNOWN" || sv_type == "NEUTRAL") {
        return;
    }
    
    if (start >= end) {
        throw std::runtime_error("ERROR: Invalid SV at position " + std::to_string(start) + "-" + std::to_string(end));
    }

    // Insert the SV call in sorted order
    SVCall sv_call{start, end, sv_type, alt_allele, data_type, genotype, hmm_likelihood, 1, 0};
    auto it = std::lower_bound(sv_calls.begin(), sv_calls.end(), sv_call);
    // sv_calls.insert(it, sv_call);

    // Update the SV type if the SV call already exists (if likelihood is
    // higher)
    if (it != sv_calls.end() && it->start == start && it->end == end)
    {
        if (hmm_likelihood != 0.0 && hmm_likelihood > it->hmm_likelihood)
        {
            // Update the SV call
            it->sv_type = sv_type;
            it->data_type = data_type;
            it->genotype = genotype;
            it->hmm_likelihood = hmm_likelihood;
            it->support++;  // Update support
        } else {
            it->support++;  // Update support
        }
    } else {
        sv_calls.insert(it, sv_call);  // Insert the new SV call
    }

    // printMessage("Adding SV call: " + std::to_string(start) + "-" + std::to_string(end) + " with length " + std::to_string(end - start) + " and type " + sv_type);
    // sv_calls.insert(SVCall{start, end, sv_type, alt_allele, data_type, genotype, hmm_likelihood, 1});
}

void updateSVType(std::vector<SVCall>& sv_calls, uint32_t start, uint32_t end, std::string sv_type, std::string data_type, std::string genotype, double hmm_likelihood)
{
    // Update the SV type for an existing SV call
    auto it = std::lower_bound(sv_calls.begin(), sv_calls.end(), SVCall{start, end, "", "", "", "", 0.0, 0, 0});
    if (it != sv_calls.end() && it->start == start && it->end == end)
    {
        it->sv_type = sv_type;
        it->data_type = data_type;
        it->genotype = genotype;
        it->hmm_likelihood = hmm_likelihood;
    } else {
        throw std::runtime_error("ERROR: SV call not found for update at position " + std::to_string(start) + "-" + std::to_string(end));
    }
}

uint32_t getSVCount(const std::vector<SVCall>& sv_calls)
{
    return (uint32_t) sv_calls.size();
}

void concatenateSVCalls(std::vector<SVCall> &target, const std::vector<SVCall>& source)
{
    // Efficiently concatenate two sets of SV calls
    // target.insert(source.begin(), source.end());
    target.insert(target.end(), source.begin(), source.end());
}

void mergeSVs(std::vector<SVCall>& sv_calls, std::unordered_map<uint32_t, uint32_t>& breakpoint_support)
{
    if (sv_calls.size() < 2) {
        return;
    }

    // Merge SV calls if they overlap
    int initial_size = sv_calls.size();
    std::vector<SVCall> merged_sv_calls;
    auto it = sv_calls.begin();
    SVCall current_merge = *it++;
    for (; it != sv_calls.end(); ++it) {
        SVCall& next = *it;

        // Find overlap
        if (next.start <= current_merge.end) {
            // Merge the SV calls if it is a subset
            if (next.end <= current_merge.end) {
                continue;
            }

            // Merge the SV calls based on HMM log likelihood (keep the higher
            // likelihood), 0.0 indicates no likelihood (Also update support)
            if (next.hmm_likelihood != 0.0) {
                if (next.hmm_likelihood > current_merge.hmm_likelihood) {
                    current_merge = next;  // Continue with the next call
                }

            // Merge based on support
            } else if (next.support > current_merge.support) {
                current_merge = next;  // Continue with the next call

            } else {
                // Merge based on breakpoint depth
                uint32_t next_depth = breakpoint_support[next.start] + breakpoint_support[next.end];
                uint32_t current_depth = breakpoint_support[current_merge.start] + breakpoint_support[current_merge.end];
                if (next_depth > current_depth) {
                    current_merge = next;  // Continue with the next call
                
                // Merge based on SV length
                } else if (next.end - next.start > current_merge.end - current_merge.start) {
                    current_merge = next;  // Continue with the next call
                }
            }

        } else {
            // No overlap: Save the previous SV and continue
            merged_sv_calls.emplace_back(current_merge);
            current_merge = next;
        }
    }

    // Add the last merged SV call
    // printMessage("Saving SV call: " + std::to_string(current_merge.start) + "-" + std::to_string(current_merge.end) + " with likelihood " + std::to_string(current_merge.hmm_likelihood));
    merged_sv_calls.emplace_back(current_merge);

    // Replace contents of the SV calls
    sv_calls = merged_sv_calls;
    
    int updated_size = sv_calls.size();
    std::cout << "Merged " << initial_size << " SV calls into " << updated_size << " SV calls" << std::endl;
}

void filterSVsWithLowSupport(std::vector<SVCall>& sv_calls, std::unordered_map<uint32_t, uint32_t>& breakpoint_support, int min_support)
{
    // Insert breakpoint support for each SV call, and remove SV calls with low
    // support
    int prev_size = sv_calls.size();
    for (auto& sv_call : sv_calls)
    {
        sv_call.total_support = breakpoint_support[sv_call.start] + breakpoint_support[sv_call.end];
        printMessage("SV call: " + std::to_string(sv_call.start) + "-" + std::to_string(sv_call.end) + " with support " + std::to_string(sv_call.total_support) + " and likelihood " + std::to_string(sv_call.hmm_likelihood) + " and length " + std::to_string(sv_call.end - sv_call.start));
    }

    // Remove SV calls with low support, unless they are large (> 20 kb)
    sv_calls.erase(std::remove_if(sv_calls.begin(), sv_calls.end(), [min_support](const SVCall& sv_call) {
        return (sv_call.total_support < min_support && (sv_call.end - sv_call.start) < 20000);
    }), sv_calls.end());

    int updated_size = sv_calls.size();
    printMessage("Filtered " + std::to_string(prev_size) + " SV calls to " + std::to_string(updated_size) + " SV calls with support >= " + std::to_string(min_support));
}
