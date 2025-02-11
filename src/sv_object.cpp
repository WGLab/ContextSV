#include "sv_object.h"

#include <algorithm>
#include <tuple>
#include <memory>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <numeric>

#include "dbscan.h"
#include "utils.h"

bool SVCall::operator<(const SVCall & other) const
{
	return start < other.start || (start == other.start && end < other.end);
}

void addSVCall(std::vector<SVCall>& sv_calls, SVCall& sv_call)
{
    if (sv_call.sv_type == SVType::UNKNOWN || sv_call.sv_type == SVType::NEUTRAL) {
        return;
    }

    // Check if the SV call is valid
    if (sv_call.start > sv_call.end) {
        printError("ERROR: Invalid SV call at position " + std::to_string(sv_call.start) + "-" + std::to_string(sv_call.end));
        return;
    }

    // Insert the SV call in sorted order
    auto it = std::lower_bound(sv_calls.begin(), sv_calls.end(), sv_call);
    sv_calls.insert(it, sv_call);

    // Insert the SV call in sorted order
    // SVCall sv_call{start, end, sv_type, alt_allele, data_type, genotype, hmm_likelihood, read_depth, 1, 1};
    // auto it = std::lower_bound(sv_calls.begin(), sv_calls.end(), sv_call);
    // sv_calls.insert(it, sv_call);
}

uint32_t getSVCount(const std::vector<SVCall>& sv_calls)
{
    return (uint32_t) sv_calls.size();
}

void concatenateSVCalls(std::vector<SVCall> &target, const std::vector<SVCall>& source)
{
    target.insert(target.end(), source.begin(), source.end());
}

void mergeSVs(std::vector<SVCall>& sv_calls, double epsilon, int min_pts)
{
    printMessage("Merging SVs with DBSCAN, eps=" + std::to_string(epsilon) + ", min_pts=" + std::to_string(min_pts));
    
    if (sv_calls.size() < 2) {
        return;
    }
    int initial_size = sv_calls.size();

    // Cluster SVs using DBSCAN for each SV type
    std::vector<SVCall> merged_sv_calls;

    // Cluster SVs using DBSCAN for each SV type
    DBSCAN dbscan(epsilon, min_pts);
    for ( const auto& sv_type : {
        SVType::DEL,
        SVType::DUP,
        SVType::INV,
        SVType::INS,
        SVType::BND,
        SVType::INV_DUP
    })
    {
        // Create a vector of SV calls for the current SV type and size interval
        std::vector<SVCall> sv_type_calls;
        std::copy_if(sv_calls.begin(), sv_calls.end(), std::back_inserter(sv_type_calls), [sv_type](const SVCall& sv_call) {
            return sv_call.sv_type == sv_type;
        });

        if (sv_type_calls.size() < 2) {
            continue;
        }

        dbscan.fit(sv_type_calls);
        const std::vector<int>& clusters = dbscan.getClusters();
        std::map<int, std::vector<SVCall>> cluster_map;
        for (size_t i = 0; i < clusters.size(); ++i) {
            cluster_map[clusters[i]].push_back(sv_type_calls[i]);
        }

        // Merge SVs in each cluster
        int cluster_count = 0;
        for (auto& cluster : cluster_map) {
            int cluster_id = cluster.first;
            std::vector<SVCall>& cluster_sv_calls = cluster.second;
            // if (cluster_id < 0) {
            //     continue;  // Skip noise and unclassified points
            // } else {
            if (true) {
                // Check if any SV has a non-zero likelihood
                bool has_nonzero_likelihood = false;
                if (cluster_sv_calls.size() > 0) {
                    for (const auto& sv_call : cluster_sv_calls) {

                        // Check if any SV has a non-zero likelihood
                        if (sv_call.hmm_likelihood != 0.0) {
                            has_nonzero_likelihood = true;
                            break;
                        }
                    }
                }
                
                SVCall merged_sv_call = cluster_sv_calls[0];
                if (has_nonzero_likelihood) {
                    // These are detected from split reads, choose the one with
                    // the highest non-zero likelihood
                    std::sort(cluster_sv_calls.begin(), cluster_sv_calls.end(), [](const SVCall& a, const SVCall& b) {
                        return a.hmm_likelihood > b.hmm_likelihood;
                    });

                    // Obtain the highest non-zero likelihood
                    auto it = std::find_if(cluster_sv_calls.begin(), cluster_sv_calls.end(), [](const SVCall& sv_call) {
                        return sv_call.hmm_likelihood != 0.0;
                    });
                    merged_sv_call = *it;

                } else {
                    // Use the median length SV
                    std::sort(cluster_sv_calls.begin(), cluster_sv_calls.end(), [](const SVCall& a, const SVCall& b) {
                        return (a.end - a.start) < (b.end - b.start);
                    });
                    int median_index = cluster_sv_calls.size() / 2;
                    merged_sv_call = cluster_sv_calls[median_index];
                }

                if (cluster_id < 0) {
                    merged_sv_call.cluster_size = cluster_id;
                } else {
                    merged_sv_call.cluster_size = (int) cluster_sv_calls.size();
                }
                merged_sv_calls.push_back(merged_sv_call);
                cluster_count++;
            }
        }
        printMessage("Completed DBSCAN with epsilon " + std::to_string(epsilon) + " for " + std::to_string(cluster_count) + " clusters of " + getSVTypeString(sv_type));
    }
    sv_calls = std::move(merged_sv_calls); // Replace with filtered list

    int updated_size = sv_calls.size();
    printMessage("Merged " + std::to_string(initial_size) + " SV calls into " + std::to_string(updated_size) + " SV calls");
}

void mergeSVSubsets(std::vector<SVCall> &sv_calls)
{
    // Sort the SV calls by start position
    int initial_size = sv_calls.size();
    std::sort(sv_calls.begin(), sv_calls.end(), [](const SVCall& a, const SVCall& b) {
        return a.start < b.start;
    });

    // Remove SVs that are subsets of other SVs
    std::vector<SVCall> filtered_sv_calls;
    // Since the input SV calls are sorted by start position, we can iterate
    // through them in order and only keep the SVs that are not subsets of
    // others
    for (const auto& sv_call : sv_calls) {
        // Check if the current SV call is a subset of any previously added
        // SV call
        bool is_subset = false;
        for (const auto& filtered_sv_call : filtered_sv_calls) {
            if (sv_call.start >= filtered_sv_call.start && sv_call.end <= filtered_sv_call.end) {
                is_subset = true;
                break;
            }
        }
        // If it's not a subset, add it to the filtered list
        if (!is_subset) {
            filtered_sv_calls.push_back(sv_call);
        }
    }
    sv_calls = std::move(filtered_sv_calls); // Replace with filtered list
    int updated_size = sv_calls.size();
    printMessage("Filtered SV calls to remove subsets, from " + std::to_string(initial_size) + " to " + std::to_string(updated_size));
}

void filterSVsWithLowSupport(std::vector<SVCall> &sv_calls, int min_support)
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
