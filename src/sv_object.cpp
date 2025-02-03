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
    sv_calls.insert(it, sv_call);

    // // Determine if the SV call already exists
    // if (it != sv_calls.end() && it->start == start && it->end == end)
    // {
    //     it->support += 1;  // Update the read support

    //     // Update SV type if likelihood is higher
    //     if (hmm_likelihood != 0.0 && hmm_likelihood > it->hmm_likelihood)
    //     {
    //         // Update the SV call
    //         it->sv_type = sv_type;
    //         it->data_type = data_type;
    //         it->genotype = genotype;
    //         it->hmm_likelihood = hmm_likelihood;
    //         it->qual = qual;
    //     }
    // } else {
    //     sv_calls.insert(it, sv_call);  // Insert the new SV call
    // }
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
    if (sv_calls.size() < 2) {
        return;
    }
    int initial_size = sv_calls.size();

    // Cluster SVs using DBSCAN for each SV type
    std::vector<SVCall> merged_sv_calls;

    // Create a set of size intervals and corresponding DBSCAN epsilons
    std::map<std::pair<int, int>, double> size_to_eps;
    // size_to_eps[{0, 1000}] = 200;
    // size_to_eps[{1000, 5000}] = 500;
    // size_to_eps[{5000, 10000}] = 3000;
    // size_to_eps[{10000, 50000}] = 4000;
    // size_to_eps[{50000, 100000}] = 5000;
    // size_to_eps[{100000, 500000}] = 10000;
    // size_to_eps[{500000, 1000000}] = 20000;

    // Small SVs
    // size_to_eps[{50, 200}] = 50;

    // // Medium SVs
    // size_to_eps[{200, 1000}] = 200;

    // // Large SVs
    // size_to_eps[{1000, 10000}] = 1000;

    // // Very large SVs
    // size_to_eps[{10000, 100000}] = 10000;

    // // Extreme SVs
    // size_to_eps[{100000, 1000000}] = 20000;

    // std::vector<double> epsilons = {50, 200, 1000, 10000, 100000, 500000};
    // std::vector<double> epsilons = {0.2};

    // double epsilon = size_interval.second;
    // // Calculate epsilon as 20% of the largest size in the interval
    // double epsilon = 0.1 * size_interval.first.second;
    // printMessage("Clustering SVs with size " + std::to_string(size_interval.first.first) + "-" + std::to_string(size_interval.first.second) + " with epsilon " + std::to_string(epsilon));
    // int min_pts = 2;
    // int min_pts = 2;
    DBSCAN dbscan(epsilon, min_pts);
    // DBSCAN dbscan(size_interval.second, 2);
    
    for ( const auto& sv_type : {
        SVType::DEL,
        SVType::DUP,
        SVType::INV,
        SVType::INS,
        SVType::BND,
        SVType::INV_DUP
    })
    {
        // DBSCAN dbscan(1000, 2);

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
            if (cluster_id < 0) {
                continue;  // Skip noise and unclassified points
            } else {
                // Use the median length SV
                std::sort(cluster_sv_calls.begin(), cluster_sv_calls.end(), [](const SVCall& a, const SVCall& b) {
                    return (a.end - a.start) < (b.end - b.start);
                });
                int median_index = cluster_sv_calls.size() / 2;
                SVCall median_sv_call = cluster_sv_calls[median_index];
                median_sv_call.cluster_size = (int) cluster_sv_calls.size();
                merged_sv_calls.push_back(median_sv_call);
                cluster_count++;
            }
        }
        printMessage("Completed DBSCAN with epsilon " + std::to_string(epsilon) + " for " + std::to_string(cluster_count) + " clusters of " + getSVTypeString(sv_type));
    }
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
