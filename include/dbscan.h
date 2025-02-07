#ifndef DBSCAN_H
#define DBSCAN_H

#include <vector>
#include <utility>
#include <cmath>
#include <algorithm>

#include "sv_object.h"

class DBSCAN {
    public:
        DBSCAN(double epsilon, int minPts) : epsilon(epsilon), minPts(minPts) {}

        // Fit the DBSCAN algorithm to SV calls
        void fit(const std::vector<SVCall>& sv_calls);

        const std::vector<int>& getClusters() const;

    private:
        double epsilon;
        int minPts;
        std::vector<int> clusters;

        // Expand the cluster for a given SV call
        bool expandCluster(const std::vector<SVCall>& sv_calls, size_t pointIdx, int clusterId);

        // Find the region query for a given SV call
        std::vector<size_t> regionQuery(const std::vector<SVCall>& sv_calls, size_t pointIdx) const;

        // Calculate the distance between two SV calls
        double distance(const SVCall& a, const SVCall& b) const;
};

#endif // DBSCAN_H
