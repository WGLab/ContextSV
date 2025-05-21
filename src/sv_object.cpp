#include "sv_object.h"

#include <algorithm>
#include <tuple>
#include <memory>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <numeric>
#include <fstream>
#include <map>

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

void mergeSVs(std::vector<SVCall>& sv_calls, double epsilon, int min_pts, bool keep_noise, const std::string& json_filepath)
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
        std::vector<SVCall> merged_sv_type_calls;

        // Create a vector of SV calls for the current SV type and size interval
        std::vector<SVCall> sv_type_calls;
        std::copy_if(sv_calls.begin(), sv_calls.end(), std::back_inserter(sv_type_calls), [sv_type](const SVCall& sv_call) {
            return sv_call.sv_type == sv_type;
        });

        if (sv_type_calls.size() < 2) {
            // Add all unclustered points to the merged list
            for (const auto& sv_call : sv_type_calls) {
                SVCall noise_sv_call = sv_call;
                // merged_sv_calls.push_back(noise_sv_call);
                merged_sv_type_calls.push_back(noise_sv_call);
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

        // Save clusters to JSON if requested
        if (!json_filepath.empty()) {
            // Create the directory if it doesn't exist
            std::string dir = json_filepath.substr(0, json_filepath.find_last_of('/'));
            if (!fileExists(dir)) {
                std::string command = "mkdir -p " + dir;
                system(command.c_str());
            }
            // Save the clusters to a JSON file
            // Prepend the SV type before the extension
            // Remove the file extension from the JSON filename
            std::string json_filename_no_ext = json_filepath.substr(0, json_filepath.find_last_of('.'));
            std::string json_filename = json_filename_no_ext + "_" + getSVTypeString(sv_type) + ".json";
            // std::string json_filename = json_filepath + "/clusters_" + getSVTypeString(sv_type) + ".json";
            saveClustersToJSON(json_filename, cluster_map);
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
                    // merged_sv_calls.push_back(noise_sv_call);
                    merged_sv_type_calls.push_back(noise_sv_call);
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
                    // merged_sv_calls.push_back(merged_sv_call);
                    merged_sv_type_calls.push_back(merged_sv_call);

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
                    // merged_sv_calls.push_back(merged_sv_call);
                    merged_sv_type_calls.push_back(merged_sv_call);
                }
                cluster_count++;
            }
        }
        printMessage("Merged " + std::to_string(cluster_count) + " clusters of " + getSVTypeString(sv_type) + ", found " + std::to_string(merged_sv_type_calls.size()) + " merged SV calls");

        // Merge overlapping SVs by cluster size
        std::sort(merged_sv_type_calls.begin(), merged_sv_type_calls.end(), [](const SVCall& a, const SVCall& b) {
            return a.start < b.start || (a.start == b.start && a.end < b.end);
        });
        std::vector<SVCall> merged_sv_calls_final;
        for (size_t i = 0; i < merged_sv_type_calls.size(); i++) {
            SVCall& sv_call = merged_sv_type_calls[i];

            // Merge cluster sizes if they overlap
            if (i > 0 && sv_call.start <= merged_sv_type_calls[i - 1].end) {
                // Keep the larger cluster size
                if (sv_call.cluster_size > merged_sv_type_calls[i - 1].cluster_size) {
                    merged_sv_calls_final.push_back(sv_call);
                }
            } else {
                merged_sv_calls_final.push_back(sv_call);
            }
        }
        printMessage("Merged " + std::to_string(merged_sv_type_calls.size()) + " overlapping SV calls into " + std::to_string(merged_sv_calls_final.size()) + " merged SV calls");
        
        // Insert merged SV calls into the final list
        merged_sv_calls.insert(merged_sv_calls.end(), merged_sv_calls_final.begin(), merged_sv_calls_final.end());

        // printMessage("Completed DBSCAN with epsilon " + std::to_string(epsilon) + " for " + std::to_string(cluster_count) + " clusters of " + getSVTypeString(sv_type) + " SVs");
    }
    sv_calls = std::move(merged_sv_calls); // Replace with filtered list
    int updated_size = sv_calls.size();
    printMessage("Merged " + std::to_string(initial_size) + " SV calls into " + std::to_string(updated_size) + " SV calls");

    // // Merge overlapping SVs by cluster size
    // std::sort(sv_calls.begin(), sv_calls.end(), [](const SVCall& a, const SVCall& b) {
    //     return a.start < b.start || (a.start == b.start && a.end < b.end);
    // });
    // std::vector<SVCall> merged_sv_calls_final;
    // for (size_t i = 0; i < sv_calls.size(); i++) {
    //     SVCall& sv_call = sv_calls[i];

    //     // Merge cluster sizes if they overlap
    //     if (i > 0 && sv_call.start <= sv_calls[i - 1].end) {
    //         // Keep the larger cluster size
    //         if (sv_call.cluster_size > sv_calls[i - 1].cluster_size) {
    //             sv_calls[i - 1] = sv_call;
    //         }
    //     } else {
    //         merged_sv_calls_final.push_back(sv_call);
    //     }
    // }
    // sv_calls = std::move(merged_sv_calls_final); // Replace with filtered list
    // int final_size = sv_calls.size();
    // printMessage("Merged " + std::to_string(updated_size) + " overlapping SV calls into " + std::to_string(final_size) + " SV calls");
}

void saveClustersToJSON(const std::string &filename, const std::map<int, std::vector<SVCall>> &clusters)
{
    // Check if the filename is empty
    if (filename.empty()) {
        printError("ERROR: Filename is empty");
        return;
    }

    // Remove the file if it already exists
    if (fileExists(filename)) {
        std::remove(filename.c_str());
    }

    // Open the JSON file for writing
    std::ofstream json_file(filename);

    if (!json_file.is_open()) {
        printError("ERROR: Unable to open JSON file for writing: " + filename);
        return;
    }
    json_file << "{\n";
    json_file << "  \"clusters\": [\n";
    size_t count = 0;
    // for (size_t i = 0; i < clusters.size(); ++i) {
    for (const auto& [cluster_id, cluster] : clusters) {
        if (cluster_id < 0) {
            continue; // Skip noise points
        }

        // const auto& cluster = clusters.at(i);
        // const auto& cluster = sv_list;
        json_file << "    {\n";
        json_file << "      \"cluster_id\": " << cluster_id << ",\n";
        json_file << "      \"cluster_size\": " << cluster.size() << ",\n";
        json_file << "      \"sv_calls\": [\n";
        for (size_t j = 0; j < cluster.size(); ++j) {
            const auto& sv_call = cluster[j];
            json_file << "        {\n";
            json_file << "          \"start\": " << sv_call.start << ",\n";
            json_file << "          \"end\": " << sv_call.end << "\n";
            // json_file << "          \"sv_type\": \"" << getSVTypeString(sv_call.sv_type) << "\",\n";
            // json_file << "          \"alt_allele\": \"" << sv_call.alt_allele << "\",\n";
            // json_file << "          \"genotype\": \"" << getGenotypeString(sv_call.genotype) << "\",\n";
            // json_file << "          \"hmm_likelihood\": " << sv_call.hmm_likelihood << "\n";
            json_file << "        }" << (j < cluster.size() - 1 ? "," : "") << "\n";
        }
        json_file << "      ]\n";
        // json_file << "    }" << (i < clusters.size() - 1 ? "," : "") << "\n";
        count++;
        if (count < clusters.size() - 1) {
            json_file << "    }," << "\n";
        } else {
            json_file << "    }\n";
            printMessage("JSON found last cluster: " + std::to_string(cluster_id));
        }
    }
    json_file << "  ]\n";
    json_file << "}\n";
    json_file.close();
    printMessage("Saved clusters to JSON file: " + filename);
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
    }
    int merge_count = initial_size - combined_sv_calls.size();
    sv_calls = std::move(combined_sv_calls); // Replace with filtered list
    if (merge_count > 0) {
        printMessage("Merged " + std::to_string(merge_count) + " SV candidates with identical start and end positions");
    }
}
