#include "dbscan1d.h"

#include <vector>
#include <cmath>
#include <algorithm>
#include <map>

void DBSCAN1D::fit(const std::vector<int>& points) {
    int clusterId = 0;
    clusters.assign(points.size(), -1); // -1 means unclassified

    for (size_t i = 0; i < points.size(); ++i) {
        if (clusters[i] == -1) { // if point is not yet classified
            if (expandCluster(points, i, clusterId)) {
                ++clusterId;
            }
        }
    }
}

const std::vector<int>& DBSCAN1D::getClusters() const {
    return clusters;
}

bool DBSCAN1D::expandCluster(const std::vector<int>& points, size_t pointIdx, int clusterId) {
    std::vector<size_t> seeds = regionQuery(points, pointIdx);
    if (static_cast<int>(seeds.size()) < minPts) {
        clusters[pointIdx] = -2; // mark as noise
        return false;
    }

    for (size_t seedIdx : seeds) {
        clusters[seedIdx] = clusterId;
    }

    seeds.erase(std::remove(seeds.begin(), seeds.end(), pointIdx), seeds.end());

    while (!seeds.empty()) {
        size_t currentPoint = seeds.back();
        seeds.pop_back();

        std::vector<size_t> result = regionQuery(points, currentPoint);
        if (static_cast<int>(result.size()) >= minPts) {
            for (size_t resultPoint : result) {
                if (clusters[resultPoint] == -1 || clusters[resultPoint] == -2) {
                    if (clusters[resultPoint] == -1) {
                        seeds.push_back(resultPoint);
                    }
                    clusters[resultPoint] = clusterId;
                }
            }
        }
    }

    return true;
}

std::vector<size_t> DBSCAN1D::regionQuery(const std::vector<int>& points, size_t pointIdx) const {
    std::vector<size_t> neighbors;
    for (size_t i = 0; i < points.size(); ++i) {
        if (distance(points[pointIdx], points[i]) <= epsilon) {
            neighbors.push_back(i);
        }
    }
    return neighbors;
}

double DBSCAN1D::distance(int point1, int point2) const {
    return std::abs(point1 - point2);
}

std::vector<int> DBSCAN1D::getLargestCluster(const std::vector<int> &points)
{
    std::vector<int> clusters = getClusters();
    std::map<int, std::vector<int>> cluster_map;
    for (size_t i = 0; i < clusters.size(); ++i) {
        cluster_map[clusters[i]].push_back(points[i]);
    }

    int largest_cluster_id = -1;
    size_t largest_size = 0;
    for (const auto &entry : cluster_map) {
        if (entry.first >= 0 && entry.second.size() > largest_size) {
            largest_size = entry.second.size();
            largest_cluster_id = entry.first;
        }
    }

    return cluster_map[largest_cluster_id];
}
