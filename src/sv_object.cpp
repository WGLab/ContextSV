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
    // Check if the SV call is valid
    if (sv_call.start > sv_call.end) {
        printError("ERROR: Invalid SV call at position " + std::to_string(sv_call.start) + "-" + std::to_string(sv_call.end) + " from data type " + getSVAlignmentTypeString(sv_call.aln_type));
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

void mergeSVs(std::vector<SVCall>& sv_calls, double epsilon, int min_pts, bool keep_noise)
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
    })
    {
        // Create a vector of SV calls for the current SV type and size interval
        std::vector<SVCall> sv_type_calls;
        std::copy_if(sv_calls.begin(), sv_calls.end(), std::back_inserter(sv_type_calls), [sv_type](const SVCall& sv_call) {
            return sv_call.sv_type == sv_type;
        });

        if (sv_type_calls.size() < 2) {
            // Add all unclustered points to the merged list
            for (const auto& sv_call : sv_type_calls) {
                SVCall noise_sv_call = sv_call;
                merged_sv_calls.push_back(noise_sv_call);
            }
            continue;
        }

        dbscan.fit(sv_type_calls);

        // Create a map of cluster IDs to SV calls
        const std::vector<int>& clusters = dbscan.getClusters();
        std::map<int, std::vector<SVCall>> cluster_map;  // Cluster ID to SV calls
        if (sv_type == SVType::INS) {
            // Add only non-CIGARCLIP SVs to the cluster map
            for (size_t i = 0; i < clusters.size(); ++i) {
                // if (sv_type_calls[i].data_type != SVDataType::CIGARCLIP) {
                // Use the SVEvidenceFlags to check for CIGARCLIP
                if (!sv_type_calls[i].aln_type.test(static_cast<size_t>(SVDataType::CIGARCLIP))) {
                    cluster_map[clusters[i]].push_back(sv_type_calls[i]);
                }
            }
        } else {
            for (size_t i = 0; i < clusters.size(); ++i) {
                cluster_map[clusters[i]].push_back(sv_type_calls[i]);
            }
        }

        // Merge SVs in each cluster
        int cluster_count = 0;
        for (auto& cluster : cluster_map) {
            int cluster_id = cluster.first;
            std::vector<SVCall>& cluster_sv_calls = cluster.second;

            // Continue if fewer than 2 SV calls in the cluster (due to CIGARCLIP filter)
            if (cluster_sv_calls.size() < 2) {
                continue;
            }

            // Add unmerged SV calls
            if (cluster_id < 0 && keep_noise) {

                // Add all unclustered points to the merged list
                for (const auto& sv_call : cluster_sv_calls) {
                    SVCall noise_sv_call = sv_call;
                    merged_sv_calls.push_back(noise_sv_call);
                }

            // Merge clustered SV calls
            } else {

                // ----------------------------
                // HMM-BASED MERGING
                // ----------------------------
                
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
                    // // These are detected from split reads, choose the one with
                    // // the highest non-zero likelihood normalized by the length of the SV
                    // std::sort(cluster_sv_calls.begin(), cluster_sv_calls.end(), [](const SVCall& a, const SVCall& b) {
                    //     return (a.hmm_likelihood / (double)(a.end - a.start + 1)) > (b.hmm_likelihood / (double)(b.end - b.start + 1));
                    // });

                    // // Obtain the highest non-zero likelihood
                    // auto it = std::find_if(cluster_sv_calls.begin(), cluster_sv_calls.end(), [](const SVCall& sv_call) {
                    //     return sv_call.hmm_likelihood != 0.0;
                    // });

                    // Choose the SV with the highest cluster size of all SVs
                    // with non-zero likelihood (if equal, choose the larger SV)
                    std::sort(cluster_sv_calls.begin(), cluster_sv_calls.end(), [](const SVCall& a, const SVCall& b) {
                        return a.cluster_size > b.cluster_size || (a.cluster_size == b.cluster_size && a.end - a.start > b.end - b.start);
                    });
                    // std::sort(cluster_sv_calls.begin(), cluster_sv_calls.end(), [](const SVCall& a, const SVCall& b) {
                    //     return a.cluster_size > b.cluster_size || (a.cluster_size == b.cluster_size && a.hmm_likelihood > b.hmm_likelihood);
                    // });
                    auto it = std::find_if(cluster_sv_calls.begin(), cluster_sv_calls.end(), [](const SVCall& sv_call) {
                        return sv_call.hmm_likelihood != 0.0;
                    });

                    // Add SV call
                    merged_sv_call = *it;
                    merged_sv_calls.push_back(merged_sv_call);

                // ----------------------------
                // CIGAR-BASED MERGING
                // ----------------------------

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

                    // Add SV call
                    merged_sv_call.cluster_size = (int) cluster_sv_calls.size();
                    merged_sv_calls.push_back(merged_sv_call);
                }
                cluster_count++;
            }
        }
        printMessage("Completed DBSCAN with epsilon " + std::to_string(epsilon) + " for " + std::to_string(cluster_count) + " clusters of " + getSVTypeString(sv_type) + " SVs");
    }
    sv_calls = std::move(merged_sv_calls); // Replace with filtered list
    int updated_size = sv_calls.size();
    printMessage("Merged " + std::to_string(initial_size) + " SV calls into " + std::to_string(updated_size) + " SV calls");
}

void mergeDuplicateSVs(std::vector<SVCall> &sv_calls)
{
    int initial_size = sv_calls.size();
    std::vector<SVCall> combined_sv_calls;

    // Sort first by start position, then by SV type
    std::sort(sv_calls.begin(), sv_calls.end(), [](const SVCall& a, const SVCall& b) {
        return std::tie(a.start, a.sv_type) < std::tie(b.start, b.sv_type);
    });
    for (size_t i = 0; i < sv_calls.size(); i++) {
        SVCall& sv_call = sv_calls[i];

        // Merge cluster sizes if start and end positions are the same
        if (i > 0 && sv_call.start == sv_calls[i - 1].start && sv_call.end == sv_calls[i - 1].end) {
            // Combine the cluster sizes
            sv_call.cluster_size += sv_calls[i - 1].cluster_size;
            combined_sv_calls.back() = sv_call;
        } else {
            combined_sv_calls.push_back(sv_call);
        }
        // SVCall& sv_call = sv_calls[i];
        // // For SVs at the same start position with the same SV type, keep the one
        // // with the highest likelihood
        // if (i > 0 && sv_call.start == sv_calls[i - 1].start && ((sv_call.sv_type == sv_calls[i - 1].sv_type) || sv_call.sv_type == SVType::UNKNOWN || sv_calls[i - 1].sv_type == SVType::UNKNOWN)) {
        //     // Keep the SV call with a non-zero likelihood
        //     // The HMM prediction is more reliable than the split read prediction
        //     if (sv_call.hmm_likelihood != 0.0 && sv_calls[i - 1].hmm_likelihood == 0.0) {
        //         // Combine the cluster sizes
        //         sv_call.cluster_size += sv_calls[i - 1].cluster_size;
        //         combined_sv_calls.back() = sv_call;
        //     }

        //     // If the likelihoods are equal, keep the one with the larger cluster size
        //     // This is to ensure that the SV call with more supporting reads is
        //     // kept
        //     else if (sv_call.hmm_likelihood == sv_calls[i - 1].hmm_likelihood && sv_call.cluster_size >= sv_calls[i - 1].cluster_size) {
        //         // Combine the cluster sizes
        //         sv_call.cluster_size += sv_calls[i - 1].cluster_size;
        //         combined_sv_calls.back() = sv_call;
        //     }
        // } else {
        //     combined_sv_calls.push_back(sv_call);
        // }
    }
    int merge_count = initial_size - combined_sv_calls.size();
    sv_calls = std::move(combined_sv_calls); // Replace with filtered list
    if (merge_count > 0) {
        printMessage("Merged " + std::to_string(merge_count) + " SV candidates with identical start and end positions");
    }
}
