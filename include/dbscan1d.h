#ifndef DBSCAN1D_H
#define DBSCAN1D_H


#include <vector>
#include <utility>
#include <cmath>
#include <algorithm>


class DBSCAN1D {
    public:
        DBSCAN1D(double epsilon, int minPts) : epsilon(epsilon), minPts(minPts) {}

        void fit(const std::vector<int>& points);

        const std::vector<int>& getClusters() const;

        std::vector<int> getLargestCluster(const std::vector<int> &points);

    private:
        double epsilon;
        int minPts;
        std::vector<int> clusters;

        bool expandCluster(const std::vector<int>& points, size_t pointIdx, int clusterId);

        std::vector<size_t> regionQuery(const std::vector<int>& points, size_t pointIdx) const;

        double distance(int a, int b) const;

};

#endif  // DBSCAN1D_H
