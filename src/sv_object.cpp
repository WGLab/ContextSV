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
    // if (sv_call.sv_type == SVType::UNKNOWN || sv_call.sv_type == SVType::NEUTRAL) {
    //     return;
    // }

    // Check if the SV call is valid
    if (sv_call.start > sv_call.end) {
        printError("ERROR: Invalid SV call at position " + std::to_string(sv_call.start) + "-" + std::to_string(sv_call.end));
        return;
    }

    // Insert the SV call in sorted order
    auto it = std::lower_bound(sv_calls.begin(), sv_calls.end(), sv_call);
    sv_calls.insert(it, sv_call);
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

    // Cluster SVs using DBSCAN for each SV type
    int initial_size = sv_calls.size();
    std::vector<SVCall> merged_sv_calls;
    DBSCAN dbscan(epsilon, min_pts);
    for ( const auto& sv_type : {
        SVType::DEL,
        SVType::DUP,
        SVType::INV,
        SVType::INS,
        SVType::BND,
        SVType::INV_DUP,
        SVType::INV_DEL,
    })
    {
        // [TEST] Skip if not insertions
        // if (sv_type != SVType::INS) {
        //     continue;
        // }

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


            // [TEST] If insertions, and if any SV has length between 9400 and
            // 9500, print all SV coordinates in the cluster
            bool print_all = false;
            // if (sv_type == SVType::INS) {
            //     for (const auto& sv_call : cluster_sv_calls) {
            //         // printMessage("[TEST] SV call " + std::to_string(sv_call.start) + "-" + std::to_string(sv_call.end) + ", length=" + std::to_string((sv_call.end - sv_call.start) + 1));
            //         // if (sv_call.end - sv_call.start >= 9400 && sv_call.end -
            //         // sv_call.start <= 9500) {
            //         // if (sv_call.end - sv_call.start >= 15100 && sv_call.end -
            //         // sv_call.start <= 15200) {
            //         // if (sv_call.end - sv_call.start >= 11200 && sv_call.end -
            //         // sv_call.start <= 11300) {
            //         // if (sv_call.end - sv_call.start >= 16800 && sv_call.end -
            //         // sv_call.start <= 17000) {
            //         // if (sv_call.end - sv_call.start >= 11300 && sv_call.end -
            //         // sv_call.start <= 11400) {
            //         // if (sv_call.end - sv_call.start >= 13100 && sv_call.end -
            //         // sv_call.start <= 13200) {
            //         if (sv_call.end - sv_call.start >= 28200 && sv_call.end - sv_call.start <= 28300) {
            //             print_all = true;
            //             break;
            //         }
            //     }
            // }
            if (print_all) {
                printMessage("[TEST] Cluster " + std::to_string(cluster_id) + " has " + std::to_string(cluster_sv_calls.size()) + " SVs:");
                for (const auto& sv_call : cluster_sv_calls) {
                    printMessage("  " + std::to_string(sv_call.start) + "-" + std::to_string(sv_call.end) + ", length=" + std::to_string((sv_call.end - sv_call.start) + 1));
                }
            }

            if (cluster_id < 0) {
                // Add all noise points to the merged list if >10 kb
                // for (const auto& sv_call : cluster_sv_calls) {
                //     if ((sv_call.end - sv_call.start)+1 >= 10000) {
                //         SVCall noise_sv_call = sv_call;
                //         noise_sv_call.cluster_size = cluster_id;
                //         merged_sv_calls.push_back(noise_sv_call);
                //         printMessage("[TEST] Adding noise SV " + std::to_string(sv_call.start) + "-" + std::to_string(sv_call.end) + ", length=" + std::to_string((sv_call.end - sv_call.start) + 1));
                //     }
                // }
                continue;  // Skip noise and unclassified points
            } else {
            // if (true) {
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

                    // [TEST]
                    if (print_all) {
                        printMessage("[TEST] Merging cluster " + std::to_string(cluster_id) + " with highest likelihood SV " + std::to_string(merged_sv_call.start) + "-" + std::to_string(merged_sv_call.end) + ", length=" + std::to_string((merged_sv_call.end - merged_sv_call.start) + 1));
                    }

                } else {
                    // Use the median length SV of the top 10% of the cluster
                    // (shorter reads are often noise)
                    std::sort(cluster_sv_calls.begin(), cluster_sv_calls.end(), [](const SVCall& a, const SVCall& b) {
                        return (a.end - a.start) > (b.end - b.start);
                    });

                    // Get the top 10% of the cluster
                    size_t top_10_percent = std::max(1, (int) (cluster_sv_calls.size() * 0.1));
                    std::vector<SVCall> top_10(cluster_sv_calls.begin(), cluster_sv_calls.begin() + top_10_percent);

                    // Get the median SV for the top 10% of the cluster
                    size_t median_index = top_10.size() / 2;
                    merged_sv_call = top_10[median_index];

                    // // Get the starting index of the top 10% of the cluster
                    // // (Cluster is sorted by descending length)
                    // size_t start_index = std::max(0, (int) (cluster_sv_calls.size() * 0.9));

                    // // Get the top 10% of the cluster
                    // std::vector<SVCall> top_half(cluster_sv_calls.begin() + start_index, cluster_sv_calls.end());

                    // // Get the median SV for the top 50% of the cluster
                    // size_t median_index = top_half.size() / 2;
                    // merged_sv_call = top_half[median_index];
                    // int median_index = cluster_sv_calls.size() / 2;
                    // merged_sv_call = cluster_sv_calls[median_index];

                    // [TEST]
                    if (print_all) {
                        printMessage("[TEST] Merging cluster " + std::to_string(cluster_id) + " with median SV " + std::to_string(merged_sv_call.start) + "-" + std::to_string(merged_sv_call.end) + ", length=" + std::to_string((merged_sv_call.end - merged_sv_call.start) + 1));
                    }
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
