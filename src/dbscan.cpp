#include "dbscan.h"

#include <vector>
#include <utility>
#include <cmath>
#include <algorithm>


void DBSCAN::fit(const std::vector<SVCall>& sv_calls) {
    int clusterId = 0;
    // clusters.assign(points.size(), -1); // -1 means unclassified
    clusters.assign(sv_calls.size(), -1); // -1 means unclassified

    // for (size_t i = 0; i < points.size(); ++i) {
    for (size_t i = 0; i < sv_calls.size(); ++i) {
        if (clusters[i] == -1) { // if point is not yet classified
            // if (expandCluster(points, i, clusterId)) {
            if (expandCluster(sv_calls, i, clusterId)) {
                ++clusterId;
            }
        }
    }
}

const std::vector<int>& DBSCAN::getClusters() const {
    return clusters;
}

// bool DBSCAN::expandCluster(const std::vector<std::pair<double, double>>&
// points, size_t pointIdx, int clusterId) {
bool DBSCAN::expandCluster(const std::vector<SVCall>& sv_calls, size_t pointIdx, int clusterId) {
    std::vector<size_t> seeds = regionQuery(sv_calls, pointIdx);
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

        std::vector<size_t> result = regionQuery(sv_calls, currentPoint);
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

std::vector<size_t> DBSCAN::regionQuery(const std::vector<SVCall>& sv_calls, size_t pointIdx) const {
    std::vector<size_t> neighbors;
    for (size_t i = 0; i < sv_calls.size(); ++i) {
        if (distance(sv_calls[pointIdx], sv_calls[i]) <= epsilon) {
            neighbors.push_back(i);
        }
    }
    return neighbors;
}

double DBSCAN::distance(const SVCall& point1, const SVCall& point2) const {
    // return std::sqrt(std::pow(point1.first - point2.first, 2) +
    // std::pow(point1.second - point2.second, 2));
    // return std::sqrt(std::pow(static_cast<double>(point1.start) - static_cast<double>(point2.start), 2) +
    // std::pow(static_cast<double>(point1.end) -
    // static_cast<double>(point2.end), 2));
    
    // Calculate reciprocal overlap-based distance
    // https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02840-6
    // https://link.springer.com/article/10.1186/gb-2009-10-10-r119
    int overlap = std::max(0, std::min(static_cast<int>(point1.end), static_cast<int>(point2.end)) - std::max(static_cast<int>(point1.start), static_cast<int>(point2.start)));
    int length1 = static_cast<int>(point1.end - point1.start);
    int length2 = static_cast<int>(point2.end - point2.start);

    // Minimum reciprocal overlap
    double distance = 1.0 - std::min(static_cast<double>(overlap) / static_cast<double>(length1), static_cast<double>(overlap) / static_cast<double>(length2));
    // double distance = 1.0 - static_cast<double>(overlap) / std::min(length1, length2);
    return distance;  // 0.0 means identical, 1.0 means no overlap
}
