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

void addSVCall(std::vector<SVCall>& sv_calls, uint32_t start, uint32_t end, std::string sv_type, const std::string& alt_allele, std::string data_type, std::string genotype, double hmm_likelihood, int read_depth)
{
    // Ignore unknown SV types
    if (sv_type == "UNKNOWN" || sv_type == "NEUTRAL") {
        return;
    }
    
    if (start > end) {
        printError("ERROR: Invalid SV at position " + std::to_string(start) + "-" + std::to_string(end));
        return;
    }

    // Insert the SV call in sorted order
    SVCall sv_call{start, end, sv_type, alt_allele, data_type, genotype, hmm_likelihood, read_depth, 1};
    auto it = std::lower_bound(sv_calls.begin(), sv_calls.end(), sv_call);

    // Update the SV type if the SV call already exists (if likelihood is
    // higher)
    if (it != sv_calls.end() && it->start == start && it->end == end)
    {
        it->support += 1;  // Update the read support
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
        printError("ERROR: SV call not found for update at position " + std::to_string(start) + "-" + std::to_string(end));
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
    // int initial_size = sv_calls.size();
    
    // Merge any SV calls that have >90% reciprocal overlap
    std::vector<SVCall> merged_sv_calls;
    SVCall current_merge = sv_calls[0];
    for (size_t i = 1; i < sv_calls.size(); i++) {
        SVCall& next = sv_calls[i];
        // Check for overlap
        if (next.start <= current_merge.end) {
            //XprintMessage("Comparing SV " + std::to_string(current_merge.start) + "-" + std::to_string(current_merge.end) + " (support " + std::to_string(current_merge.support) + ", length " + std::to_string(current_merge.end - current_merge.start) + ") with " + std::to_string(next.start) + "-" + std::to_string(next.end) + " (support " + std::to_string(next.support) + ", length " + std::to_string(next.end - next.start) + ")");
            
            // if (current_merge.start <= next.end && next.start <= current_merge.end) {
            // Calculate reciprocal overlap
            uint32_t overlap = std::max(0, (int)std::min(current_merge.end, next.end) - (int)std::max(current_merge.start, next.start));
            uint32_t union_length = std::max(current_merge.end, next.end) - std::min(current_merge.start, next.start);
            double overlap_fraction = static_cast<double>(overlap) / union_length;
            //XprintMessage("Overlap fraction: " + std::to_string(overlap_fraction));

            // Merge if reciprocal overlap is >90%
            if (overlap_fraction > 0.90) {
                //XprintMessage("Merging SV calls with overlap " + std::to_string(overlap_fraction));
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
                }
            } else {
            	// Continue with the larger length
				uint32_t current_length = current_merge.end - current_merge.start;
				uint32_t next_length = next.end - next.start;
				if (next_length > current_length) {  // And support meets threshold
					current_merge = next;
				}
            }
        } else {
            merged_sv_calls.push_back(current_merge);
            current_merge = next;
        }
    }

    merged_sv_calls.push_back(current_merge);  // Add the last SV call
    sv_calls = merged_sv_calls;  // Update the SV calls

    // int updated_size = sv_calls.size();
    // std::cout << "Merged " << initial_size << " SV calls into " << updated_size << " SV calls" << std::endl;
}

void filterSVsWithLowSupport(std::vector<SVCall>& sv_calls, int min_support)
{
    // Filter SV calls with low read support
    sv_calls.erase(std::remove_if(sv_calls.begin(), sv_calls.end(), [min_support](const SVCall& sv_call) {
        return sv_call.support < min_support;
    }), sv_calls.end());
}

void filterSVsWithLowSupport(std::vector<SVCall> &sv_calls, int min_support, const std::string &data_type)
{
    // Filter SV calls with low read depth only for the specified data type, keeping the rest
    sv_calls.erase(std::remove_if(sv_calls.begin(), sv_calls.end(), [min_support, data_type](const SVCall& sv_call) {
        return sv_call.support < min_support && sv_call.data_type == data_type;
    }), sv_calls.end());
}
