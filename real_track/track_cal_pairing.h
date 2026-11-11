#pragma once

#include <vector>
#include <utility>

#include "track.h"
#include "cluster.h"

using namespace std;

vector<pair<track*, const cluster*>> matchTracksToClusters(
    vector<track>& realTracks,
    const vector<cluster>& clusters
);

void calculateAndSetMass(
    track& realTrack,
    const cluster& calorimeterCluster
){
    realTrack.setMass()
}