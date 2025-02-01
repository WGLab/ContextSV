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

        void fit(const std::vector<SVCall>& sv_calls);

        const std::vector<int>& getClusters() const;

    private:
        double epsilon;
        int minPts;
        std::vector<int> clusters;

        bool expandCluster(const std::vector<SVCall>& sv_calls, size_t pointIdx, int clusterId);

        std::vector<size_t> regionQuery(const std::vector<SVCall>& sv_calls, size_t pointIdx) const;

        double distance(const SVCall& a, const SVCall& b) const;
};

#endif // DBSCAN_H
