#ifndef KALMAN_H
#define KALMAN_H

#include "track.h"
#include "settings.h"   // ZOrder

#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <limits>

using namespace std;

// =====================================================
// ===== 2D KALMAN STATE (one projection: ZX or ZY) ====
// =====================================================
/**
 * State of a track in one projection:
 *   coord = measured coordinate (x for ZX, y for ZY)
 *   slope = derivative d(coord)/dz
 *
 * Covariance terms:
 *   varCoord       = Var(coord)
 *   varSlope       = Var(slope)
 *   covCoordSlope  = Cov(coord, slope)
 */
struct KalmanState2D {
    double coord;
    double slope;
    double varCoord;
    double varSlope;
    double covCoordSlope;
};

// =====================================================
// ===== TRACK FINDING CONFIGURATION ===================
// =====================================================
struct CKFConfig {
    int minHits = 3;
    double chi2Cut = 9.0;
    double processNoiseSlope = 0.0;
    double chi2AmbiguityCut = 1.0;
    bool greedy = true;
    bool generous = false;
    ZOrder zOrder = ZOrder::DecreasingZ;

    void setGreedy(bool value) { greedy = value; }

    void setGenerous(bool value) {
        generous = value;
        if (value) greedy = false;
    }
};

// =====================================================
// ===== INITIALISE KALMAN STATE FROM TWO HITS =========
// =====================================================
/**
 * Initialises the state using two hits (seed).
 * The order of the two hits follows the configured ZOrder.
 */
inline KalmanState2D initKalmanState(
    const Hit* firstHit,
    const Hit* secondHit,
    PlaneType axis)
{
    const double z1 = firstHit->getZ();
    const double z2 = secondHit->getZ();
    const double coord1 = firstHit->getCoord(axis);
    const double coord2 = secondHit->getCoord(axis);
    const double sigma1 = firstHit->getSigma();
    const double sigma2 = secondHit->getSigma();
    const double dz = z2 - z1;

    KalmanState2D state;
    state.slope = (coord2 - coord1) / dz;
    state.coord = coord2;
    state.varCoord = sigma2 * sigma2;
    state.varSlope = (sigma1 * sigma1 + sigma2 * sigma2) / (dz * dz);
    state.covCoordSlope = (sigma2 * sigma2) / dz;
    return state;
}

// =====================================================
// ===== PREDICT: EXTRAPOLATE TO NEXT PLANE =============
// =====================================================
inline void predictKalmanState(
    KalmanState2D& state,
    double dz,
    double processNoiseSlope = 0.0)
{
    const double coordPred = state.coord + state.slope * dz;
    const double slopePred = state.slope;
    const double varCoordPred = state.varCoord
                              + dz * dz * state.varSlope
                              + 2.0 * dz * state.covCoordSlope;
    const double covPred = state.covCoordSlope + dz * state.varSlope;
    const double varSlopePred = state.varSlope + processNoiseSlope;

    state.coord = coordPred;
    state.slope = slopePred;
    state.varCoord = varCoordPred;
    state.covCoordSlope = covPred;
    state.varSlope = varSlopePred;
}

// =====================================================
// ===== CHI² FOR A CANDIDATE HIT (GATING) ==============
// =====================================================
inline double computeChi2(
    const KalmanState2D& state,
    double measuredCoord,
    double measurementSigma)
{
    const double residualVariance = state.varCoord
                                  + measurementSigma * measurementSigma;
    const double residual = measuredCoord - state.coord;
    return (residual * residual) / residualVariance;
}

// =====================================================
// ===== UPDATE: CORRECT STATE WITH THE CHOSEN HIT ======
// =====================================================
inline void updateKalmanState(
    KalmanState2D& state,
    double measuredCoord,
    double measurementSigma)
{
    const double residualVariance = state.varCoord
                                  + measurementSigma * measurementSigma;
    const double kalmanGainCoord = state.varCoord / residualVariance;
    const double kalmanGainSlope = state.covCoordSlope / residualVariance;
    const double residual = measuredCoord - state.coord;

    const double coordUpd = state.coord + kalmanGainCoord * residual;
    const double slopeUpd = state.slope + kalmanGainSlope * residual;
    const double varCoordUpd = (1.0 - kalmanGainCoord) * state.varCoord;
    const double covUpd = (1.0 - kalmanGainCoord) * state.covCoordSlope;
    const double varSlopeUpd = state.varSlope
                             - (state.covCoordSlope * state.covCoordSlope)
                               / residualVariance;

    state.coord = coordUpd;
    state.slope = slopeUpd;
    state.varCoord = varCoordUpd;
    state.covCoordSlope = covUpd;
    state.varSlope = varSlopeUpd;
}

// =====================================================
// ===== GROUP HITS BY PLANE ===========================
// =====================================================
/**
 * Groups hits by their Z-coordinate and keeps only hits measuring axis.
 * The resulting planes are ordered according to zOrder.
 */
inline vector<vector<const Hit*> > groupHitsByPlane(
    const vector<const Hit*>& hits,
    PlaneType axis,
    ZOrder zOrder,
    bool alreadyGrouped = false)
{
    (void)alreadyGrouped; // reserved for a future optimisation
    map<double, vector<const Hit*> > hitsGroupedByZ;

    for (const Hit* h : hits) {
        if (h && h->hasCoord(axis))
            hitsGroupedByZ[h->getZ()].push_back(h);
    }

    vector<vector<const Hit*> > hitsByPlane;
    hitsByPlane.reserve(hitsGroupedByZ.size());

    for (const auto& zAndHits : hitsGroupedByZ)
        hitsByPlane.push_back(zAndHits.second);

    if (zOrder == ZOrder::DecreasingZ)
        reverse(hitsByPlane.begin(), hitsByPlane.end());

    return hitsByPlane;
}

// Overload for the vector<Hit*> used by the event loop.
inline vector<vector<const Hit*> > groupHitsByPlane(
    const vector<Hit*>& hits,
    PlaneType axis,
    ZOrder zOrder,
    bool alreadyGrouped = false)
{
    vector<const Hit*> constHits(hits.begin(), hits.end());
    return groupHitsByPlane(constHits, axis, zOrder, alreadyGrouped);
}

// =====================================================
// ===== BUILD A SINGLE TRACK FROM A SEED (GREEDY) ======
// =====================================================
/**
 * Runs the CKF from a seed through all remaining planes.
 * On every plane the hit with the smallest acceptable chi2 is selected.
 * A candidate is discarded if a plane has no acceptable hit.
 */
inline vector<const Hit*> buildTrackFromSeed(
    const vector<vector<const Hit*> >& hitsByPlane,
    PlaneType axis,
    const Hit* firstHit,
    const Hit* secondHit,
    const CKFConfig cfg,
    int& ambiguitySkipped)
{
    vector<const Hit*> trackHits = { firstHit, secondHit };
    KalmanState2D state = initKalmanState(firstHit, secondHit, axis);
    double previousZ = secondHit->getZ();

    const int nPlanes = static_cast<int>(hitsByPlane.size());
    for (int planeIndex = 2; planeIndex < nPlanes; ++planeIndex) {
        const vector<const Hit*>& hitsOnThisPlane = hitsByPlane[planeIndex];
        if (hitsOnThisPlane.empty())
            return {};

        const double currentZ = hitsOnThisPlane.front()->getZ();
        const double dz = currentZ - previousZ;
        predictKalmanState(state, dz, cfg.processNoiseSlope);

        const Hit* bestHit = nullptr;
        double bestChi2 = numeric_limits<double>::infinity();
        double secondBestChi2 = numeric_limits<double>::infinity();

        for (const Hit* candidateHit : hitsOnThisPlane) {
            const double candidateChi2 = computeChi2(
                state, candidateHit->getCoord(axis), candidateHit->getSigma());

            if (candidateChi2 < bestChi2) {
                secondBestChi2 = bestChi2;
                bestChi2 = candidateChi2;
                bestHit = candidateHit;
            } else if (candidateChi2 < secondBestChi2) {
                secondBestChi2 = candidateChi2;
            }
        }

        if (bestHit == nullptr || bestChi2 >= cfg.chi2Cut)
            return {};

        if (secondBestChi2 < numeric_limits<double>::infinity()) {
            const double diff = secondBestChi2 - bestChi2;
            if (diff < cfg.chi2AmbiguityCut) {
                ++ambiguitySkipped;
                return {};
            }
        }

        updateKalmanState(state, bestHit->getCoord(axis), bestHit->getSigma());
        trackHits.push_back(bestHit);
        previousZ = currentZ;
    }

    if (static_cast<int>(trackHits.size()) < cfg.minHits)
        return {};

    return trackHits;
}

// =====================================================
// ===== RECURSIVE TRACK BUILDING (GENEROUS MODE) =======
// =====================================================
inline void buildTracksRecursive(
    const vector<vector<const Hit*> >& hitsByPlane,
    PlaneType axis,
    const CKFConfig& cfg,
    int planeIndex,
    double previousZ,
    KalmanState2D state,
    vector<const Hit*>& currentTrack,
    vector<vector<const Hit*> >& result,
    int& ambiguitySkipped)
{
    if (ambiguitySkipped > 0)
        return;

    if (planeIndex == static_cast<int>(hitsByPlane.size())) {
        if (static_cast<int>(currentTrack.size()) >= cfg.minHits)
            result.push_back(currentTrack);
        return;
    }

    const vector<const Hit*>& hitsOnThisPlane = hitsByPlane[planeIndex];
    if (hitsOnThisPlane.empty())
        return;

    const double currentZ = hitsOnThisPlane.front()->getZ();
    const double dz = currentZ - previousZ;
    predictKalmanState(state, dz, cfg.processNoiseSlope);

    vector<pair<const Hit*, double> > candidates;
    double bestChi2 = numeric_limits<double>::infinity();
    double secondBestChi2 = numeric_limits<double>::infinity();

    for (const Hit* h : hitsOnThisPlane) {
        const double candidateChi2 = computeChi2(
            state, h->getCoord(axis), h->getSigma());
        if (candidateChi2 >= cfg.chi2Cut)
            continue;

        candidates.push_back(make_pair(h, candidateChi2));
        if (candidateChi2 < bestChi2) {
            secondBestChi2 = bestChi2;
            bestChi2 = candidateChi2;
        } else if (candidateChi2 < secondBestChi2) {
            secondBestChi2 = candidateChi2;
        }
    }

    if (candidates.empty())
        return;

    if (secondBestChi2 < numeric_limits<double>::infinity() &&
        secondBestChi2 - bestChi2 < cfg.chi2AmbiguityCut) {
        ++ambiguitySkipped;
        return;
    }

    for (const auto& candidate : candidates) {
        KalmanState2D nextState = state;
        updateKalmanState(nextState,
                          candidate.first->getCoord(axis),
                          candidate.first->getSigma());
        currentTrack.push_back(candidate.first);

        buildTracksRecursive(
            hitsByPlane, axis, cfg,
            planeIndex + 1,
            currentZ,
            nextState,
            currentTrack,
            result,
            ambiguitySkipped);

        currentTrack.pop_back();
    }
}

// Build all possible tracks from one seed (generous mode).
inline vector<vector<const Hit*> > buildTracksFromSeedGenerous(
    const vector<vector<const Hit*> >& hitsByPlane,
    PlaneType axis,
    const Hit* firstHit,
    const Hit* secondHit,
    const CKFConfig cfg,
    int& ambiguitySkipped)
{
    vector<vector<const Hit*> > result;
    vector<const Hit*> currentTrack = { firstHit, secondHit };
    KalmanState2D state = initKalmanState(firstHit, secondHit, axis);

    buildTracksRecursive(
        hitsByPlane, axis, cfg,
        2,
        secondHit->getZ(),
        state,
        currentTrack,
        result,
        ambiguitySkipped);

    return result;
}

// =====================================================
// ===== MAIN ITERATIVE TRACK FINDER (CKF) =============
// =====================================================
/**
 * Finds projection candidates with the Combinatorial Kalman Filter.
 * Each returned candidate is a vector of hit pointers in one projection.
 */
inline vector<vector<const Hit*> > findTracksCKF(
    vector<vector<const Hit*> > hitsByPlane,
    PlaneType axis,
    int& ambiguitySkipped,
    const CKFConfig cfg)
{
    vector<vector<const Hit*> > foundTracks;
    if (hitsByPlane.size() < 2)
        return foundTracks;

    // ---- Generous mode ----
    if (cfg.generous) {
        for (const Hit* firstHit : hitsByPlane[0]) {
            for (const Hit* secondHit : hitsByPlane[1]) {
                vector<vector<const Hit*> > candidates =
                    buildTracksFromSeedGenerous(
                        hitsByPlane, axis, firstHit, secondHit,
                        cfg, ambiguitySkipped);

                if (ambiguitySkipped > 0)
                    return {};

                foundTracks.insert(foundTracks.end(),
                                   candidates.begin(), candidates.end());
            }
        }
        return foundTracks;
    }

    // ---- Greedy mode ----
    while (true) {
        vector<const Hit*> foundTrack;

        for (const Hit* firstHit : hitsByPlane[0]) {
            for (const Hit* secondHit : hitsByPlane[1]) {
                vector<const Hit*> candidate = buildTrackFromSeed(
                    hitsByPlane, axis, firstHit, secondHit,
                    cfg, ambiguitySkipped);

                if (ambiguitySkipped > 0)
                    return {};

                if (!candidate.empty()) {
                    foundTrack = candidate;
                    break;
                }
            }
            if (!foundTrack.empty())
                break;
        }

        if (foundTrack.empty())
            break;

        foundTracks.push_back(foundTrack);

        if (cfg.greedy) {
            for (size_t planeIndex = 0; planeIndex < foundTrack.size(); ++planeIndex) {
                vector<const Hit*>& plane = hitsByPlane[planeIndex];
                plane.erase(remove(plane.begin(), plane.end(),
                                   foundTrack[planeIndex]), plane.end());
            }
        }
    }

    return foundTracks;
}

#endif // KALMAN_H
