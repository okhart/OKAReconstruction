#ifndef CLUSTER_BUILDER_H
#define CLUSTER_BUILDER_H

#include <vector>
#include "Hit.h"
#include "Cluster.h"

class ClusterBuilder {
public:
    ClusterBuilder();
    ~ClusterBuilder() = default;

    void SetHits(const std::vector<CalHit>& hits);  //Load all hits of a detector

    void BuildClusters();   //Make separate clusters

    const std::vector<Cluster>& GetClusters() const;    //return clusters

    void Clear();   //Make ready for the next event

private:
    std::vector<CalHit> hits;
    std::vector<Cluster> clusters;
};

#endif