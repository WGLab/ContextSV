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

        std::vector<SVCall> cluster;
        cluster.push_back(sv_calls[i]);
        merged[i] = true;

        // Use 10% of the length of the first SV as the threshold
        uint32_t sv_a_window = (uint32_t) std::ceil((double) (sv_calls[i].end - sv_calls[i].start + 1) * 0.1);

        // Find SVs that have start or end positions within 10% of each other's length
        for (size_t j = i + 1; j < sv_calls.size(); j++) {
            if (merged[j]) continue;

            // Check if the SVs are within 10% of the largest SV's length
            uint32_t sv_b_window = (uint32_t) std::ceil((double) (sv_calls[j].end - sv_calls[j].start + 1) * 0.1);
            uint32_t sv_window = std::max(sv_a_window, sv_b_window);
            bool start_within_window = std::abs((int) sv_calls[j].start - (int) sv_calls[i].start) <= (int) sv_window;
            bool end_within_window = std::abs((int) sv_calls[j].end - (int) sv_calls[i].end) <= (int) sv_window;
            if (start_within_window && end_within_window) {
                cluster.push_back(sv_calls[j]);
                merged[j] = true;
            }
        }

        // Remove clusters with single SVs that have low support
        if (cluster.size() < 2 && cluster[0].support < 2) {
            continue;
        }

        std::vector<SVCall> filtered_cluster = cluster;

        // If any SV length equals 2039, print all the SV calls in the cluster
        bool found_2039 = false;
        for (const auto& sv : filtered_cluster) {
            if (sv.end - sv.start == 2039) {
                printMessage("[TEST] Found SV with length 2039 at " + std::to_string(sv.start) + "-" + std::to_string(sv.end) + " (SUP=" + std::to_string(sv.support) + ")");
                found_2039 = true;
            }
        }
        if (found_2039) {
            std::cout << "[TEST] Cluster of SVs with size " << filtered_cluster.size() << ":" << std::endl;
            for (const auto& sv : filtered_cluster) {
                printMessage("SV: " + std::to_string(sv.start) + "-" + std::to_string(sv.end) + " (SUP=" + std::to_string(sv.support) + ", LEN=" + std::to_string(sv.end - sv.start) + ")");
            }
        }

        // Find the median-length SV in the cluster and use it as the merged SV
        // Sort the cluster by length
        std::sort(filtered_cluster.begin(), filtered_cluster.end(), [](const SVCall& a, const SVCall& b) {
            return (a.end - a.start) < (b.end - b.start);
        });

        // Get the median SV
        size_t median_index = filtered_cluster.size() / 2;
        SVCall median_sv = filtered_cluster[median_index];

        median_sv.cluster_size = (int) cluster.size();
        
        // Add the merged SV to the list
        merged_sv_calls.push_back(median_sv);
    }
    sv_calls = std::move(merged_sv_calls); // Replace with filtered list

    // Print SVs that have length 2039
    for (const auto& sv_call : sv_calls) {
        if (sv_call.end - sv_call.start == 2039) {
            printMessage("[TEST] Found merged SV with length 2039 at " + std::to_string(sv_call.start) + "-" + std::to_string(sv_call.end) + " (SUP=" + std::to_string(sv_call.support) + ")");
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
