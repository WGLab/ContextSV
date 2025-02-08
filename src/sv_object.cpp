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

void addSVCall(std::vector<SVCall>& sv_calls, uint32_t start, uint32_t end, SVType sv_type, const std::string& alt_allele, std::string data_type, std::string genotype, double hmm_likelihood, int read_depth)
{
    if (sv_type == SVType::UNKNOWN || sv_type == SVType::NEUTRAL) {
        return;
    }
    
    if (start > end) {
        printError("ERROR: Invalid SV at position " + std::to_string(start) + "-" + std::to_string(end));
        return;
    }

    // Insert the SV call in sorted order
    SVCall sv_call{start, end, sv_type, alt_allele, data_type, genotype, hmm_likelihood, read_depth, 1, 1};
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
        // Create a DBSCAN object for the current SV type
        // epsilon = 0.45;
        // min_pts = 15;
        // if (sv_type == SVType::DEL) {
        //     epsilon = 0.45;
        //     min_pts = 16;
        // } else {
        //     // epsilon = 0.65;
        //     // min_pts = 15;
        //     // epsilon = 0.45;
        //     // min_pts = 16;
        //     // epsilon = 0.45;
        //     // min_pts = 2;
        //     // epsilon = 0.45;
        //     // min_pts = 15;
        // }
        // DBSCAN dbscan(epsilon, min_pts);

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
                // Use the highest HMM likelihood normalized by SV size as the
                // representative SV (if any non-zero likelihoods exist)
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

                // [TEST] Check if any SV has a length greater than 600kb
                bool found_large_sv = false;
                for (const auto& sv_call : cluster_sv_calls) {
                    if (sv_call.end - sv_call.start > 600000) {
                        found_large_sv = true;
                        break;
                    }
                }
                if (found_large_sv) {
                    printMessage("Found large SV with length greater than 600kb");
                    printMessage("Found " + std::to_string(cluster_sv_calls.size()) + " SVs in cluster " + std::to_string(cluster_id) + " of type " + getSVTypeString(sv_type) + ", with epsilon=" + std::to_string(epsilon) + ", min_pts=" + std::to_string(min_pts));
                }
                
                SVCall merged_sv_call = cluster_sv_calls[0];
                if (has_nonzero_likelihood) {
                    // Use the highest HMM likelihood normalized by SV size as the
                    // representative SV
                    // std::vector<double> likelihoods;
                    // Default very low log-likelihood for zero likelihoods
                    std::vector<double> likelihoods(cluster_sv_calls.size(), -std::numeric_limits<double>::infinity());
                    // for (const auto& sv_call : cluster_sv_calls) {
                    int i = 0;
                    for (const auto& sv_call : cluster_sv_calls) {
                        if (sv_call.hmm_likelihood != 0.0) {
                            uint32_t sv_size = (uint32_t) (sv_call.end - sv_call.start);
                            if (sv_size > 0) {
                                likelihoods[i] = sv_call.hmm_likelihood / sv_size;
                                // likelihoods.push_back(sv_call.hmm_likelihood / sv_size);
                            }
                        }

                        // Print the SV length, likelihood, and normalized
                        // likelihood
                        if (found_large_sv) {
                            printMessage("Start: " + std::to_string(sv_call.start) + ", end: " + std::to_string(sv_call.end) + ", likelihood: " + std::to_string(sv_call.hmm_likelihood) + ", normalized likelihood: " + std::to_string(likelihoods[i]) + ", length: " + std::to_string(sv_call.end - sv_call.start));
                            // printMessage("SV length: " + std::to_string(sv_call.end - sv_call.start) + ", likelihood: " + std::to_string(sv_call.hmm_likelihood) + ", normalized likelihood: " + std::to_string(likelihoods[i]) + ", start: " + std::to_string(sv_call.start) + ", end: " + std::to_string(sv_call.end));
                        }
                        i++;
                    }
                    
                    // Find the index of the maximum element in the likelihoods
                    // vector
                    auto max_likelihood_it = std::max_element(likelihoods.begin(), likelihoods.end());
                    int max_likelihood_index = std::distance(likelihoods.begin(), max_likelihood_it);
                    merged_sv_call = cluster_sv_calls[max_likelihood_index];
                    printMessage("Merged SV with highest normalized likelihood: " + std::to_string(merged_sv_call.start) + "-" + std::to_string(merged_sv_call.end) + ", likelihood: " + std::to_string(merged_sv_call.hmm_likelihood) + ", normalized likelihood: " + std::to_string(merged_sv_call.hmm_likelihood / (merged_sv_call.end - merged_sv_call.start)) + ", size: " + std::to_string(merged_sv_call.end - merged_sv_call.start));

                } else {
                    // Use the median length SV
                    std::sort(cluster_sv_calls.begin(), cluster_sv_calls.end(), [](const SVCall& a, const SVCall& b) {
                        return (a.end - a.start) < (b.end - b.start);
                    });
                    int median_index = cluster_sv_calls.size() / 2;
                    merged_sv_call = cluster_sv_calls[median_index];
                    printMessage("Merged SV with median length: " + std::to_string(merged_sv_call.start) + "-" + std::to_string(merged_sv_call.end) + ", likelihood: " + std::to_string(merged_sv_call.hmm_likelihood) + ", size: " + std::to_string(merged_sv_call.end - merged_sv_call.start));
                
                if (cluster_id < 0) {
                    merged_sv_call.cluster_size = cluster_id;
                } else {
                    merged_sv_call.cluster_size = (int) cluster_sv_calls.size();
                }
                // merged_sv_call.cluster_size = (int) cluster_sv_calls.size();
                merged_sv_calls.push_back(merged_sv_call);
                cluster_count++;
                }
            }
        }
        printMessage("Completed DBSCAN with epsilon " + std::to_string(epsilon) + " for " + std::to_string(cluster_count) + " clusters of " + getSVTypeString(sv_type));
    }

    printMessage("[TEST] Merged " + std::to_string(initial_size) + " SV calls into " + std::to_string(merged_sv_calls.size()) + " SV calls");
    sv_calls = std::move(merged_sv_calls); // Replace with filtered list

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
