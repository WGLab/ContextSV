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

void addSVCall(std::vector<SVCall>& sv_calls, uint32_t start, uint32_t end, std::string sv_type, std::string alt_allele, std::string data_type, std::string genotype, double hmm_likelihood, int read_depth)
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
    SVCall sv_call{start, end, sv_type, alt_allele, data_type, genotype, hmm_likelihood, read_depth, 1};
    auto it = std::lower_bound(sv_calls.begin(), sv_calls.end(), sv_call);

    // Update the SV type if the SV call already exists (if likelihood is
    // higher)
    if (it != sv_calls.end() && it->start == start && it->end == end)
    {
        it->support += 1;  // Update the read support
        // printMessage("Updating SV call with length " + std::to_string(end - start) + " and type " + sv_type + " and support " + std::to_string(it->support));
        if (hmm_likelihood != 0.0 && hmm_likelihood > it->hmm_likelihood)
        {
            // Update the SV call
            it->sv_type = sv_type;
            it->data_type = data_type;
            it->genotype = genotype;
            it->hmm_likelihood = hmm_likelihood;
        }
    } else {
        sv_calls.insert(it, sv_call);  // Insert the new SV call
    }
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
    target.insert(target.end(), source.begin(), source.end());
}

void mergeSVs(std::vector<SVCall>& sv_calls)
{
    if (sv_calls.size() < 2) {
        return;
    }

    // Merge SV calls if they overlap
    int initial_size = sv_calls.size();
    
    // Merge any SV calls that have >90% reciprocal overlap
    std::vector<SVCall> merged_sv_calls;
    SVCall current_merge = sv_calls[0];
    for (size_t i = 1; i < sv_calls.size(); i++) {
        SVCall& next = sv_calls[i];
        // Check for overlap
        if (next.start <= current_merge.end) {
            // printMessage("Comparing SV " + std::to_string(current_merge.start) + "-" + std::to_string(current_merge.end) + " (support " + std::to_string(current_merge.support) + ", length " + std::to_string(current_merge.end - current_merge.start) + ") with " + std::to_string(next.start) + "-" + std::to_string(next.end) + " (support " + std::to_string(next.support) + ", length " + std::to_string(next.end - next.start) + ")");
            // if (current_merge.start <= next.end && next.start <= current_merge.end) {
            // Calculate reciprocal overlap
            uint32_t overlap = std::max(0, (int)std::min(current_merge.end, next.end) - (int)std::max(current_merge.start, next.start));
            uint32_t union_length = std::max(current_merge.end, next.end) - std::min(current_merge.start, next.start);
            double overlap_fraction = static_cast<double>(overlap) / union_length;
            // printMessage("Overlap fraction: " + std::to_string(overlap_fraction));

            // Merge if reciprocal overlap is >90%
            if (overlap_fraction > 0.90) {
                // printMessage("Merging SV calls with overlap " + std::to_string(overlap_fraction));
                // Keep the SV call with the higher read support
                if (next.support > current_merge.support) {
                    current_merge = next;
                } else if (next.support == current_merge.support) {
                    // Keep the SV call with the higher likelihood
                    if (next.hmm_likelihood != 0.0 && current_merge.hmm_likelihood != 0.0 && next.hmm_likelihood > current_merge.hmm_likelihood) {
                        current_merge = next;
                    } else if (next.hmm_likelihood == current_merge.hmm_likelihood) {
                        // Keep the SV call with the higher read depth
                        if (next.read_depth > current_merge.read_depth) {
                            current_merge = next;
                        }
                    }
                    // // Keep the SV call with the higher read depth
                    // if (next.read_depth > current_merge.read_depth) {
                    //     current_merge = next;
                    // } else if (next.read_depth == current_merge.read_depth) {
                    //     // Keep the SV call with the higher likelihood
                    //     if (next.hmm_likelihood > current_merge.hmm_likelihood) {
                    //         current_merge = next;
                    //     }
                    // }
                }
            } else {
                merged_sv_calls.push_back(current_merge);
                current_merge = next;
            }
        } else {
            merged_sv_calls.push_back(current_merge);
            current_merge = next;
        }
    }

    // Add the last SV call
    merged_sv_calls.push_back(current_merge);

    // Update the SV calls
    sv_calls = merged_sv_calls;
    // for (size_t i = 0; i < sv_calls.size(); i++) {
    //     SVCall& current = sv_calls[i];
    //     bool merged = false;
    //     for (size_t j = i + 1; j < sv_calls.size(); j++) {
    //         SVCall& next = sv_calls[j];
    //         if (current.start <= next.end && next.start <= current.end) {
    //             // Calculate reciprocal overlap
    //             uint32_t overlap = std::max(0, (int)std::min(current.end, next.end) - (int)std::max(current.start, next.start));
    //             uint32_t union_length = std::max(current.end, next.end) - std::min(current.start, next.start);
    //             double overlap_fraction = static_cast<double>(overlap) / union_length;

    //             // Merge if reciprocal overlap is >90%
    //             if (overlap_fraction > 0.9) {
    //                 // Keep the SV call with the higher likelihood
    //                 if (next.hmm_likelihood > current.hmm_likelihood) {
    //                     current = next;
    //                 }
    //                 merged = true;
    //             }

    //             // Remove the merged SV call
    //             sv_calls.erase(sv_calls.begin() + j);
    //             j--;

    //     }
    //     if (!merged) {
    //         merged_sv_calls.push_back(current);
    //     }
    // }
    


    // std::vector<SVCall> merged_sv_calls;
    // auto it = sv_calls.begin();
    // SVCall current_merge = *it++;
    // double log_lh_eps = 1.0;  // Log likelihood epsilon
    // for (; it != sv_calls.end(); ++it) {
    //     SVCall& next = *it;

    //     // Find overlap
    //     // printMessage("[0] Current SV call: " + std::to_string(current_merge.start) + "-" + std::to_string(current_merge.end) + " with likelihood " + std::to_string(current_merge.hmm_likelihood) + " and read depth " + std::to_string(current_merge.read_depth) + " and length " + std::to_string(current_merge.end - current_merge.start) + " and support " + std::to_string(current_merge.support));
    //     // printMessage("[0] Next SV call: " + std::to_string(next.start) + "-" + std::to_string(next.end) + " with likelihood " + std::to_string(next.hmm_likelihood) + " and read depth " + std::to_string(next.read_depth) + " and length " + std::to_string(next.end - next.start) + " and support " + std::to_string(next.support));
    //     if (next.start <= current_merge.end) {

    //         // Merge based on read support
    //         if (next.support > current_merge.support) {
    //             // Compare only if lengths are within 20% of each other
    //             uint32_t current_length = current_merge.end - current_merge.start;
    //             uint32_t next_length = next.end - next.start;
    //             double length_diff = std::abs((int)current_length - (int)next_length);
    //             double length_threshold = 0.2 * (int)current_length;
    //             if (length_diff <= length_threshold) {
    //                 current_merge = next;  // Continue with the next call
    //                 // printMessage("Keeping next SV call with support " + std::to_string(next.support));
    //             } else {
    //                 // Keep the larger SV
    //                 if (next_length > current_length) {
    //                     current_merge = next;
    //                     // printMessage("Keeping next SV call with length " + std::to_string(next_length));
    //                 }
    //             }
    //             // printMessage("Keeping next SV call with support " + std::to_string(next.support));

    //         } else if (next.support == current_merge.support) {
    //             // Merge based on existence of predictions
    //             if (next.hmm_likelihood != 0.0 && current_merge.hmm_likelihood == 0.0) {
    //                 current_merge = next;  // Continue with the next call
    //                 // printMessage("Keeping next SV call with likelihood " + std::to_string(next.hmm_likelihood));

    //             // Merge based on prediction log likelihood
    //             } else if (next.hmm_likelihood != 0.0 && current_merge.hmm_likelihood != 0.0) {
                    
    //                 // Print all SV information
    //                 // printMessage("Current SV call: " + std::to_string(current_merge.start) + "-" + std::to_string(current_merge.end) + " with likelihood " + std::to_string(current_merge.hmm_likelihood) + " and read depth " + std::to_string(current_merge.read_depth) + " and length " + std::to_string(current_merge.end - current_merge.start) + " and support " + std::to_string(current_merge.support));
    //                 // printMessage("Next SV call: " + std::to_string(next.start) + "-" + std::to_string(next.end) + " with likelihood " + std::to_string(next.hmm_likelihood) + " and read depth " + std::to_string(next.read_depth) + " and length " + std::to_string(next.end - next.start) + " and support " + std::to_string(next.support));
    //                 // printMessage("Comparing likelihoods: " + std::to_string(current_merge.hmm_likelihood) + " vs " + std::to_string(next.hmm_likelihood));

    //                 // Keep the SV call with the higher likelihood. Compare only if
    //                 // lengths are within 20% of each other
    //                 uint32_t current_length = current_merge.end - current_merge.start;
    //                 uint32_t next_length = next.end - next.start;
    //                 double length_diff = std::abs((int)current_length - (int)next_length);
    //                 double length_threshold = 0.2 * (int)current_length;
    //                 if (length_diff <= length_threshold) {
    //                     // printMessage("Length difference is within threshold: " + std::to_string(length_diff) + " <= " + std::to_string(length_threshold));

    //                     if (next.hmm_likelihood > current_merge.hmm_likelihood) {
    //                         current_merge = next;  // Continue with the next call
    //                         // printMessage("Keeping next SV call with likelihood " + std::to_string(next.hmm_likelihood));
    //                     }
                    
    //                 } else {
    //                     // Keep the larger SV
    //                     if (next_length > current_length) {
    //                         current_merge = next;
    //                         // printMessage("[2] Keeping next SV call with length " + std::to_string(next_length));
    //                     }
    //                 }
    //             }
    //         }

    //     } else {
    //         // No overlap: Save the call and continue
    //         merged_sv_calls.emplace_back(current_merge);
    //         current_merge = next;
    //     }
    // }
    // merged_sv_calls.emplace_back(current_merge);  // Save the last call
    // sv_calls = merged_sv_calls;  // Update the SV calls

    int updated_size = sv_calls.size();
    std::cout << "Merged " << initial_size << " SV calls into " << updated_size << " SV calls" << std::endl;
}

void filterSVsWithLowSupport(std::vector<SVCall>& sv_calls, int min_support)
{
    int prev_size = sv_calls.size();

    // Filter SV calls with low read support
    sv_calls.erase(std::remove_if(sv_calls.begin(), sv_calls.end(), [min_support](const SVCall& sv_call) {
        return sv_call.support < min_support;
    }), sv_calls.end());

    // // Print read depth for each SV call
    // for (const auto& sv_call : sv_calls) {
    //     std::cout << "SV call: " << sv_call.start << "-" << sv_call.end << " with depth " << sv_call.read_depth << " and length " << (sv_call.end - sv_call.start) << std::endl;
    // }

    // // Remove SV calls with low read depth
    // sv_calls.erase(std::remove_if(sv_calls.begin(), sv_calls.end(), [min_depth](const SVCall& sv_call) {
    //     return sv_call.read_depth < min_depth;
    // }), sv_calls.end());

    // int updated_size = sv_calls.size();
    // printMessage("Filtered " + std::to_string(prev_size) + " SV calls to " + std::to_string(updated_size) + " SV calls with DP >= " + std::to_string(min_depth));
}
