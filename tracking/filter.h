#ifndef FILTER_H
#define FILTER_H

#include <vector>
#include <algorithm>

#include "track.h"

using namespace std;

// Compare tracks by reconstruction quality:
// more hits first, then lower chi2.
inline bool compareTracksByQuality(const track& a, const track& b) {
    if (a.getNHits() != b.getNHits())
        return a.getNHits() > b.getNHits();
    return a.getChi2() < b.getChi2();
}

// Count hits shared by two tracks.
inline int countCommonHits(const track& t1, const track& t2) {
    int count = 0;

    for (int i = 0; i < t1.getNHits(); ++i) {
        for (int j = 0; j < t2.getNHits(); ++j) {
            if (t1.getHit(i) == t2.getHit(j)) {
                ++count;
                break;
            }
        }
    }

    return count;
}

// Remove duplicate / overlapping tracks.
// The input is not modified; the returned tracks are independent values.
inline vector<track> filterTracks(const vector<track>& allTracks) {
    vector<track> validTracks;
    validTracks.reserve(allTracks.size());

    for (const track& t : allTracks) {
        if (t.getNHits() >= 2)
            validTracks.push_back(t);
    }

    sort(validTracks.begin(), validTracks.end(), compareTracksByQuality);

    vector<track> filtered;
    filtered.reserve(validTracks.size());

    for (const track& candidate : validTracks) {
        bool hasOverlap = false;

        for (const track& accepted : filtered) {
            if (countCommonHits(candidate, accepted) >= 1) {
                hasOverlap = true;
                break;
            }
        }

        if (!hasOverlap)
            filtered.push_back(candidate);
    }

    return filtered;
}

#endif // FILTER_H
