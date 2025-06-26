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
#include "debug.h"

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

    // Set this to print cluster information for a specific SV call for debugging
    // This is useful for debugging purposes to see how the SVs are merged
    bool debug_mode = false;
    int debug_start = 10414914;  // Set to -1 to disable
    int debug_svlen_min = 15000;
    int debug_svlen_max = 16000;
    SVType debug_sv_type = SVType::INV;

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
        // Skip if not the debug SV type
        if (debug_mode && (sv_type != debug_sv_type)) {
            DEBUG_PRINT("DEBUG: Skipping SV type " + getSVTypeString(sv_type) + " for debug mode");
            continue;
        }

        DEBUG_PRINT("Merging SV type: " + getSVTypeString(sv_type) + " (epsilon=" + std::to_string(epsilon) + ", min_pts=" + std::to_string(min_pts) + ", num SVs=" + std::to_string(sv_calls.size()) + ")");
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
                merged_sv_type_calls.push_back(noise_sv_call);
            }
            continue;
        }

        dbscan.fit(sv_type_calls);

        // Create a map of cluster IDs to SV calls
        const std::vector<int>& clusters = dbscan.getClusters();
        std::map<int, std::vector<SVCall>> cluster_map;  // Cluster ID to SV calls
        for (size_t i = 0; i < clusters.size(); ++i) {
            cluster_map[clusters[i]].push_back(sv_type_calls[i]);
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
            saveClustersToJSON(json_filename, cluster_map);
        }

        // Merge SVs in each cluster
        int cluster_count = 0;
        for (auto& cluster : cluster_map) {
            int cluster_id = cluster.first;
            std::vector<SVCall>& cluster_sv_calls = cluster.second;

            // Continue unless the debug SV call is in the cluster
            // if (debug_mode && cluster_id >= 0) {
            //     if (!cluster_sv_calls.empty() &&
            //         std::any_of(cluster_sv_calls.begin(), cluster_sv_calls.end(),
            //             [debug_start, debug_sv_type, debug_svlen_min, debug_svlen_max](const SVCall& sv_call) {
            //                 const int len = std::abs(static_cast<int>(sv_call.end - sv_call.start));

            //                 const bool start_ok = (debug_start < 0 || static_cast<int>(sv_call.start) == debug_start);

            //                 const bool len_ok = (debug_svlen_min == -1 || len >= debug_svlen_min) &&
            //                                     (debug_svlen_max == -1 || len <= debug_svlen_max);

            //                 const bool type_ok = (debug_sv_type == SVType::UNKNOWN || sv_call.sv_type == debug_sv_type);

            //                 return start_ok && len_ok && type_ok;
            //             }
            //         )) {
            //         DEBUG_PRINT("DEBUG: Found SV call in noise cluster " + std::to_string(cluster_id) + " with type " + getSVTypeString(debug_sv_type));

            //     } else {
            //         continue;
            //     }
            // }

            // Continue if fewer than 2 SV calls in the cluster (due to CIGARCLIP filter)
            if (cluster_sv_calls.size() < 2) {
                continue;
            }

            // Add unmerged SV calls
            if (cluster_id < 0 && keep_noise) {

                // Add all unclustered points to the merged list
                for (const auto& sv_call : cluster_sv_calls) {
                    SVCall noise_sv_call = sv_call;
                    merged_sv_type_calls.push_back(noise_sv_call);

                    // Print the added SV calls if >10 kb and the debug SV type
                    if (debug_mode && noise_sv_call.sv_type == debug_sv_type && (noise_sv_call.end - noise_sv_call.start) > 10000) {
                        DEBUG_PRINT("DEBUG: Adding noise SV call at " + std::to_string(noise_sv_call.start) + "-" + std::to_string(noise_sv_call.end) +
                                    ", type: " + getSVTypeString(noise_sv_call.sv_type) +
                                    ", length: " + std::to_string(noise_sv_call.end - noise_sv_call.start) +
                                    ", cluster size: " + std::to_string(noise_sv_call.cluster_size) +
                                    ", likelihood: " + std::to_string(noise_sv_call.hmm_likelihood));
                    }
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
                    auto it = std::find_if(cluster_sv_calls.begin(), cluster_sv_calls.end(), [](const SVCall& sv_call) {
                        return sv_call.hmm_likelihood != 0.0;
                    });

                    // Add SV call
                    merged_sv_call = *it;
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

                    // Print the added SV calls if >10 kb and the debug SV type
                    if (debug_mode && sv_type == debug_sv_type) {
                        DEBUG_PRINT("DEBUG: Cluster " + std::to_string(cluster_id) + " with " + std::to_string(cluster_sv_calls.size()) + " SV calls (length sorted):");
                        for (const auto& sv_call : cluster_sv_calls) {
                            if ((sv_call.end - sv_call.start) > 10000) {
                                DEBUG_PRINT("DEBUG: SV call at " + std::to_string(sv_call.start) + "-" + std::to_string(sv_call.end) +
                                            ", type: " + getSVTypeString(sv_call.sv_type) +
                                            ", length: " + std::to_string(sv_call.end - sv_call.start) +
                                            ", cluster size: " + std::to_string(sv_call.cluster_size) +
                                            ", likelihood: " + std::to_string(sv_call.hmm_likelihood));
                            }
                        }
                    }

                    // Get the top % of the cluster
                    double top_pct = 0.2;
                    size_t top_pct_size = std::max(1, (int) (cluster_sv_calls.size() *  top_pct));
                    std::vector<SVCall> top_pct_calls(cluster_sv_calls.begin(), cluster_sv_calls.begin() + top_pct_size);

                    // Print the added SV calls if >10 kb and the debug SV type
                    if (debug_mode && sv_type == debug_sv_type) {
                        DEBUG_PRINT("DEBUG: Top  " + std::to_string((int)(top_pct * 100)) + "% of cluster " + std::to_string(cluster_id) + " with " +
                                    std::to_string(top_pct_calls.size()) + " SV calls (length sorted):");
                        for (const auto& sv_call : top_pct_calls) {
                            if ((sv_call.end - sv_call.start) > 10000) {
                                DEBUG_PRINT("DEBUG: SV call at " + std::to_string(sv_call.start) + "-" + std::to_string(sv_call.end) +
                                            ", type: " + getSVTypeString(sv_call.sv_type) +
                                            ", length: " + std::to_string(sv_call.end - sv_call.start) +
                                            ", cluster size: " + std::to_string(sv_call.cluster_size) +
                                            ", likelihood: " + std::to_string(sv_call.hmm_likelihood));
                            }
                        }
                    }

                    // Get the median SV for the top % of the cluster
                    size_t median_index = top_pct_calls.size() / 2;
                    merged_sv_call = top_pct_calls[median_index];

                    // Print the merged SV call
                    if (debug_mode && sv_type == debug_sv_type) {
                        DEBUG_PRINT("DEBUG: Merged SV call at " + std::to_string(merged_sv_call.start) + "-" + std::to_string(merged_sv_call.end) +
                                    ", type: " + getSVTypeString(merged_sv_call.sv_type) +
                                    ", length: " + std::to_string(merged_sv_call.end - merged_sv_call.start) +
                                    ", cluster size: " + std::to_string(merged_sv_call.cluster_size) +
                                    ", likelihood: " + std::to_string(merged_sv_call.hmm_likelihood));
                    }

                    // Add SV call
                    merged_sv_call.cluster_size = (int) cluster_sv_calls.size();
                    merged_sv_type_calls.push_back(merged_sv_call);
                }
                cluster_count++;
            }
        }
        DEBUG_PRINT("Merged " + std::to_string(cluster_count) + " clusters of " + getSVTypeString(sv_type) + ", found " + std::to_string(merged_sv_type_calls.size()) + " merged SV calls");

        // Print SV call start, end, type, and length for debugging if > 10 kb
        if (debug_mode && sv_type == debug_sv_type) {
            DEBUG_PRINT("DEBUG: Merged SV calls for " + getSVTypeString(sv_type) + ":");
            for (const auto& sv_call : merged_sv_type_calls) {
                // if ((int)sv_call.start == debug_start) {
                if ((sv_call.end - sv_call.start) > 10000) {
                    DEBUG_PRINT("DEBUG: SV call at " + std::to_string(sv_call.start) + "-" + std::to_string(sv_call.end) +
                                ", type: " + getSVTypeString(sv_call.sv_type) +
                                ", length: " + std::to_string(sv_call.end - sv_call.start) +
                                ", cluster size: " + std::to_string(sv_call.cluster_size) +
                                ", likelihood: " + std::to_string(sv_call.hmm_likelihood));
                }
            }
        }

        /*
        // Merge overlapping SVs by SV length
        std::sort(merged_sv_type_calls.begin(), merged_sv_type_calls.end(), [](const SVCall& a, const SVCall& b) {
            return a.start < b.start || (a.start == b.start && a.end < b.end);
        });
        std::vector<SVCall> merged_sv_calls_final;
        for (size_t i = 0; i < merged_sv_type_calls.size(); i++) {
            SVCall& sv_call = merged_sv_type_calls[i];

            // Merge cluster sizes if they overlap
            if (i > 0 && sv_call.start <= merged_sv_type_calls[i - 1].end) {
                // Keep the larger SV call (end - start) if they overlap
                if ((sv_call.end - sv_call.start) > (merged_sv_type_calls[i - 1].end - merged_sv_type_calls[i - 1].start)) {
                    merged_sv_type_calls[i - 1] = sv_call; // Replace the previous SV call with the current one
                }
                // Keep the larger cluster size
                // if (sv_call.cluster_size > merged_sv_type_calls[i - 1].cluster_size) {
                //     merged_sv_calls_final.push_back(sv_call);
                // }
            } else {
                merged_sv_calls_final.push_back(sv_call);
            }
        }
        DEBUG_PRINT("Merged " + std::to_string(merged_sv_type_calls.size()) + " overlapping SV calls into " + std::to_string(merged_sv_calls_final.size()) + " merged SV calls");
        
        // Print merged SV calls for debugging
        if (debug_mode) {
            DEBUG_PRINT("DEBUG: Final merged SV calls for " + getSVTypeString(sv_type) + ":");
            for (const auto& sv_call : merged_sv_calls_final) {
                // if ((int)sv_call.start == debug_start) {
                if (sv_call.sv_type == SVType::DUP) {
                    DEBUG_PRINT("DEBUG: SV call at " + std::to_string(sv_call.start) + "-" + std::to_string(sv_call.end) +
                                ", type: " + getSVTypeString(sv_call.sv_type) +
                                ", length: " + std::to_string(sv_call.end - sv_call.start) +
                                ", cluster size: " + std::to_string(sv_call.cluster_size) +
                                ", likelihood: " + std::to_string(sv_call.hmm_likelihood));
                }
            }
        }

        // Insert merged SV calls into the final list
        merged_sv_calls.insert(merged_sv_calls.end(),
        merged_sv_calls_final.begin(), merged_sv_calls_final.end());
        */
        merged_sv_calls.insert(merged_sv_calls.end(),
                               merged_sv_type_calls.begin(), merged_sv_type_calls.end());
    }
    sv_calls = std::move(merged_sv_calls); // Replace with filtered list
    int updated_size = sv_calls.size();
    printMessage("Merged " + std::to_string(initial_size) + " SV calls into " + std::to_string(updated_size) + " SV calls");
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
    for (const auto& [cluster_id, cluster] : clusters) {
        if (cluster_id < 0) {
            continue; // Skip noise points
        }

        json_file << "    {\n";
        json_file << "      \"cluster_id\": " << cluster_id << ",\n";
        json_file << "      \"cluster_size\": " << cluster.size() << ",\n";
        json_file << "      \"sv_calls\": [\n";
        for (size_t j = 0; j < cluster.size(); ++j) {
            const auto& sv_call = cluster[j];
            json_file << "        {\n";
            json_file << "          \"start\": " << sv_call.start << ",\n";
            json_file << "          \"end\": " << sv_call.end << "\n";
            json_file << "        }" << (j < cluster.size() - 1 ? "," : "") << "\n";
        }
        json_file << "      ]\n";
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
