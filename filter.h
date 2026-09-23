#ifndef FILTER_H
#define FILTER_H

#include <TFile.h>
#include <TTree.h>
#include <TVector.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <TH2F.h>
#include "track.h"

using namespace std;

// =====================================================
// ===== TRACK COMPARISON FUNCTION =====================
// =====================================================

/**
 * Compares two tracks by quality:
 * - First by number of hits (descending)
 * - Then by χ² (ascending)
 */
bool compareTracksByQuality(track* a, track* b) {
    if (!a || !b) return false;
    if (a->getNHits() != b->getNHits())
        return a->getNHits() > b->getNHits();
    return a->getChi2() < b->getChi2();
}

// =====================================================
// ===== COUNT COMMON HITS =============================
// =====================================================

/// Returns the number of hits shared between two tracks.
int countCommonHits(track* t1, track* t2) {
    if (!t1 || !t2) return 0;
    int count = 0;
    for (int i = 0; i < t1->getNHits(); i++) {
        for (int j = 0; j < t2->getNHits(); j++) {
            if (t1->getHit(i) == t2->getHit(j)) {
                count++;
                break;
            }
        }
    }
    return count;
}

// =====================================================
// ===== TRACK FILTER FUNCTION =========================
// =====================================================

/**
 * Filters overlapping tracks to remove duplicates and bad tracks.
 * Strategy:
 *   1. Remove tracks with < 2 hits.
 *   2. Sort by quality (more hits better, then lower χ²).
 *   3. Greedy selection: keep a track only if it shares < 2 hits
 *      with any already-kept track.
 *
 * @param allTracks  Input vector of tracks (will be emptied).
 * @return           Filtered vector of unique tracks.
 */
vector<track*> filterTracks(vector<track*>& allTracks) {
    if (allTracks.empty()) return {};

    // Remove invalid tracks (null or insufficient hits)
    vector<track*> validTracks;
    for (track* t : allTracks) {
        if (t && t->getNHits() >= 2) {
            validTracks.push_back(t);
        } else {
            delete t;
        }
    }
    allTracks.clear();

    if (validTracks.empty()) return {};

    // Sort by quality
    sort(validTracks.begin(), validTracks.end(), compareTracksByQuality);

    // Greedy selection
    vector<track*> filtered;
    vector<bool> used(validTracks.size(), false);

    for (size_t i = 0; i < validTracks.size(); i++) {
        if (used[i]) continue;

        bool hasOverlap = false;
        for (size_t j = 0; j < filtered.size(); j++) {
            if (countCommonHits(validTracks[i], filtered[j]) >= 1) {
                hasOverlap = true;
                break;
            }
        }

        if (!hasOverlap) {
            filtered.push_back(validTracks[i]);
            used[i] = true;
        }
    }

    // Delete unselected tracks
    for (size_t i = 0; i < validTracks.size(); i++) {
        if (!used[i]) {
            delete validTracks[i];
        }
    }

    return filtered;
}

#endif // FILTER_H