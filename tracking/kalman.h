#ifndef KALMAN_H
#define KALMAN_H

#include "track.h"
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>

using namespace std;

// =====================================================
// ===== 2D KALMAN STATE (one projection: ZX or ZY) ====
// =====================================================

/**
 * State of a track in a single projection (ZX or ZY):
 *   coord   = measured coordinate (x for ZX, y for ZY)
 *   slope   = derivative d(coord)/dz
 *
 * Covariance matrix (uncertainties):
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
    int minHits = 3;               ///< Minimum hits for a valid track
    double chi2Cut = 9.0;          ///< Maximum χ² for hit acceptance
    double processNoiseSlope = 0.0;///< Process noise added to slope (multiple scattering)
    double chi2AmbiguityCut = 1.0; ///< If (second-best - best) < this, reject as ambiguous
    bool greedy = true;            ///< Greedy: remove hits after each track
    bool generous = false;         ///< Generous: explore all branches (overrides greedy)

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
 * Initialises the Kalman state using two hits (seed).
 *
 * Derivation (error propagation):
 *   slope0 = (coord2 - coord1) / dz
 *   coord0 = coord2
 *
 *   Var(coord0)          = sigma2^2
 *   Var(slope0)          = (sigma1^2 + sigma2^2) / dz^2
 *   Cov(coord0, slope0)  = sigma2^2 / dz
 */
inline KalmanState2D initKalmanState(const hit* firstHit, const hit* secondHit, PlaneType axis) {
    double z1 = firstHit->getZ();
    double z2 = secondHit->getZ();
    double coord1 = firstHit->getCoord(axis);
    double coord2 = secondHit->getCoord(axis);
    double sigma1 = firstHit->getSigma();
    double sigma2 = secondHit->getSigma();
    double dz = z2 - z1;

    KalmanState2D state;
    state.slope = (coord2 - coord1) / dz;
    state.coord = coord2;

    state.varCoord      = sigma2 * sigma2;
    state.varSlope      = (sigma1 * sigma1 + sigma2 * sigma2) / (dz * dz);
    state.covCoordSlope = (sigma2 * sigma2) / dz;

    return state;
}

// =====================================================
// ===== PREDICT: extrapolate to next plane ============
// =====================================================

/**
 * Predicts the state at the next measurement plane.
 *
 *   coord_pred = coord + slope * dz
 *   slope_pred = slope
 *
 * Covariance propagation:
 *   Var(coord_pred) = Var(coord) + dz^2*Var(slope) + 2*dz*Cov(coord,slope)
 *   Var(slope_pred) = Var(slope) + processNoiseSlope
 *   Cov(coord_pred, slope_pred) = Cov(coord,slope) + dz*Var(slope)
 */
inline void predictKalmanState(KalmanState2D& state, double dz, double processNoiseSlope = 0.0) {
    double coordPred = state.coord + state.slope * dz;
    double slopePred = state.slope;

    double varCoordPred = state.varCoord
                         + dz * dz * state.varSlope
                         + 2.0 * dz * state.covCoordSlope;

    double covPred = state.covCoordSlope + dz * state.varSlope;
    double varSlopePred = state.varSlope + processNoiseSlope;

    state.coord = coordPred;
    state.slope = slopePred;
    state.varCoord = varCoordPred;
    state.covCoordSlope = covPred;
    state.varSlope = varSlopePred;
}

// =====================================================
// ===== CHI² FOR A CANDIDATE HIT (gating) =============
// =====================================================

/**
 * Computes the χ² contribution of a hit candidate.
 *
 *   residual = measuredCoord - coord_pred
 *   S        = Var(coord_pred) + sigma_hit^2
 *   chi2     = residual^2 / S
 */
inline double computeChi2(const KalmanState2D& state, double measuredCoord, double measurementSigma) {
    double residualVariance = state.varCoord + measurementSigma * measurementSigma;
    double residual = measuredCoord - state.coord;
    return (residual * residual) / residualVariance;
}

// =====================================================
// ===== UPDATE: correct state with the chosen hit =====
// =====================================================

/**
 * Updates the Kalman state after associating a hit.
 *
 * Kalman gain:
 *   S               = Var(coord_pred) + sigma_hit^2
 *   K_coord         = Var(coord_pred) / S
 *   K_slope         = Cov(coord_pred, slope_pred) / S
 *
 *   coord_upd = coord_pred + K_coord * residual
 *   slope_upd = slope_pred + K_slope * residual
 *
 * Covariance update:
 *   Var(coord_upd)       = (1 - K_coord) * Var(coord_pred)
 *   Cov(coord,slope)_upd = (1 - K_coord) * Cov(coord_pred,slope_pred)
 *   Var(slope_upd)       = Var(slope_pred) - Cov(coord_pred,slope_pred)^2 / S
 */
inline void updateKalmanState(KalmanState2D& state, double measuredCoord, double measurementSigma) {
    double residualVariance = state.varCoord + measurementSigma * measurementSigma;
    double kalmanGainCoord = state.varCoord / residualVariance;
    double kalmanGainSlope = state.covCoordSlope / residualVariance;
    double residual = measuredCoord - state.coord;

    double coordUpd = state.coord + kalmanGainCoord * residual;
    double slopeUpd = state.slope + kalmanGainSlope * residual;

    double varCoordUpd = (1.0 - kalmanGainCoord) * state.varCoord;
    double covUpd = (1.0 - kalmanGainCoord) * state.covCoordSlope;
    double varSlopeUpd = state.varSlope - (state.covCoordSlope * state.covCoordSlope) / residualVariance;

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
 * Groups hits by their Z-coordinate (plane), keeping only those that
 * measure the requested axis.
 *
 * @param hits           Input hits.
 * @param axis           Projection (ZX or ZY).
 * @param alreadyGrouped If true, assumes hits are already sorted by Z.
 * @return               Vector of planes, each plane is a vector of hits.
 */
inline vector<vector<const hit*>> groupHitsByPlane(
    const vector<const hit*>& hits, PlaneType axis, bool alreadyGrouped = false)
{
    (void)alreadyGrouped;  // placeholder for future optimisation

    map<double, vector<const hit*>> hitsGroupedByZ;

    for (const hit* h : hits) {
        if (h->hasCoord(axis)) {
            double z = h->getZ();
            hitsGroupedByZ[z].push_back(h);
        }
    }

    vector<vector<const hit*>> hitsByPlane;
    hitsByPlane.reserve(hitsGroupedByZ.size());
    for (const auto& zAndHits : hitsGroupedByZ) {
        hitsByPlane.push_back(zAndHits.second);
    }
    return hitsByPlane;
}

/// Overload for non-const hit pointers.
inline vector<vector<const hit*>> groupHitsByPlane(
    const vector<hit*>& hits, PlaneType axis, bool alreadyGrouped = false)
{
    vector<const hit*> constHits(hits.begin(), hits.end());
    return groupHitsByPlane(constHits, axis, alreadyGrouped);
}

// =====================================================
// ===== BUILD A SINGLE TRACK FROM A SEED (GREEDY) =====
// =====================================================

/**
 * Runs the CKF from a seed (firstHit, secondHit) through all remaining planes,
 * taking on each plane the hit with the smallest χ².
 *
 * If no acceptable hit is found on any plane, the track is discarded
 * (no missing planes allowed).
 *
 * @param hitsByPlane     Pre-grouped hits by plane.
 * @param axis            Projection.
 * @param firstHit        First hit of the seed.
 * @param secondHit       Second hit of the seed.
 * @param cfg             CKF configuration.
 * @param ambiguitySkipped Output: incremented if ambiguity is detected.
 * @return                Vector of hits forming the track, or empty if failed.
 */
inline vector<const hit*> buildTrackFromSeed(
    const vector<vector<const hit*>>& hitsByPlane,
    PlaneType axis,
    const hit* firstHit,
    const hit* secondHit,
    const CKFConfig cfg,
    int& ambiguitySkipped)
{
    vector<const hit*> trackHits = { firstHit, secondHit };

    KalmanState2D state = initKalmanState(firstHit, secondHit, axis);
    double previousZ = secondHit->getZ();

    int nPlanes = hitsByPlane.size();
    for (int planeIndex = 2; planeIndex < nPlanes; ++planeIndex) {
        const vector<const hit*>& hitsOnThisPlane = hitsByPlane[planeIndex];
        if (hitsOnThisPlane.empty()) {
            return {};
        }

        double currentZ = hitsOnThisPlane.front()->getZ();
        double dz = currentZ - previousZ;
        predictKalmanState(state, dz, cfg.processNoiseSlope);

        // Find best and second-best hits on this plane
        const hit* bestHit = nullptr;
        double bestChi2 = numeric_limits<double>::infinity();
        double secondBestChi2 = numeric_limits<double>::infinity();

        for (const hit* candidateHit : hitsOnThisPlane) {
            double chi2 = computeChi2(state, candidateHit->getCoord(axis), candidateHit->getSigma());
            if (chi2 < bestChi2) {
                secondBestChi2 = bestChi2;
                bestChi2 = chi2;
                bestHit = candidateHit;
            } else if (chi2 < secondBestChi2) {
                secondBestChi2 = chi2;
            }
        }

        if (bestHit == nullptr || bestChi2 >= cfg.chi2Cut) {
            return {};
        }

        // Check for ambiguity (two hits with similar χ²)
        if (secondBestChi2 < numeric_limits<double>::infinity()) {
            double diff = secondBestChi2 - bestChi2;
            if (diff < cfg.chi2AmbiguityCut) {
                ambiguitySkipped++;
                return {};
            }
        }

        updateKalmanState(state, bestHit->getCoord(axis), bestHit->getSigma());
        trackHits.push_back(bestHit);
        previousZ = currentZ;
    }

    if ((int)trackHits.size() < cfg.minHits) {
        return {};
    }
    return trackHits;
}

// =====================================================
// ===== RECURSIVE TRACK BUILDING (GENEROUS MODE) ======
// =====================================================

/**
 * Recursively explores all possible hit combinations (branching) for a seed.
 * Used only in generous mode.
 */
inline void buildTracksRecursive(
    const vector<vector<const hit*>>& hitsByPlane,
    PlaneType axis,
    const CKFConfig& cfg,
    int planeIndex,
    double previousZ,
    KalmanState2D state,
    vector<const hit*>& currentTrack,
    vector<vector<const hit*>>& result,
    int& ambiguitySkipped)
{
    if (ambiguitySkipped > 0)
        return;

    if (planeIndex == (int)hitsByPlane.size()) {
        if ((int)currentTrack.size() >= cfg.minHits)
            result.push_back(currentTrack);
        return;
    }

    const vector<const hit*>& hitsOnThisPlane = hitsByPlane[planeIndex];
    if (hitsOnThisPlane.empty())
        return;

    double currentZ = hitsOnThisPlane.front()->getZ();
    double dz = currentZ - previousZ;
    predictKalmanState(state, dz, cfg.processNoiseSlope);

    vector<pair<const hit*, double>> candidates;  // (hit, χ²)
    double bestChi2 = numeric_limits<double>::infinity();
    double secondBestChi2 = numeric_limits<double>::infinity();

    for (const hit* h : hitsOnThisPlane) {
        double chi2 = computeChi2(state, h->getCoord(axis), h->getSigma());
        if (chi2 >= cfg.chi2Cut)
            continue;

        candidates.push_back({h, chi2});

        if (chi2 < bestChi2) {
            secondBestChi2 = bestChi2;
            bestChi2 = chi2;
        } else if (chi2 < secondBestChi2) {
            secondBestChi2 = chi2;
        }
    }

    if (candidates.empty())
        return;

    // Ambiguity check at this plane
    if (secondBestChi2 < numeric_limits<double>::infinity() &&
        secondBestChi2 - bestChi2 < cfg.chi2AmbiguityCut)
    {
        ambiguitySkipped++;
        return;
    }

    // Branch over all viable candidates
    for (auto& candidate : candidates) {
        KalmanState2D nextState = state;
        updateKalmanState(nextState, candidate.first->getCoord(axis), candidate.first->getSigma());

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

/**
 * Wrapper for generous mode: builds all possible tracks from a seed.
 */
inline vector<vector<const hit*>> buildTracksFromSeedGenerous(
    const vector<vector<const hit*>>& hitsByPlane,
    PlaneType axis,
    const hit* firstHit,
    const hit* secondHit,
    const CKFConfig cfg,
    int& ambiguitySkipped)
{
    vector<vector<const hit*>> result;
    vector<const hit*> currentTrack = { firstHit, secondHit };

    KalmanState2D state = initKalmanState(firstHit, secondHit, axis);

    buildTracksRecursive(
        hitsByPlane, axis, cfg,
        2,                 // start from plane index 2 (third plane)
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
 * Main function to find tracks using the Combinatorial Kalman Filter.
 *
 * In each round, it tries all seeds (pairs of hits on the first two planes).
 * When a valid track is found (greedy mode), those hits are removed and the
 * process repeats. In generous mode, all possible branches are explored.
 *
 * @param hitsByPlane      Pre-grouped hits (will be modified in greedy mode).
 * @param axis             Projection (ZX or ZY).
 * @param ambiguitySkipped Output: number of tracks rejected due to ambiguity.
 * @param cfg              CKF configuration.
 * @return                 Vector of found tracks (each as a vector of hit pointers).
 */
inline vector<vector<const hit*>> findTracksCKF(
    vector<vector<const hit*>> hitsByPlane,
    PlaneType axis,
    int& ambiguitySkipped,
    const CKFConfig cfg)
{
    vector<vector<const hit*>> foundTracks;

    if (hitsByPlane.size() < 2) {
        return foundTracks;
    }

    // =====================================================
    // ===== GENEROUS MODE ================================
    // =====================================================
    if (cfg.generous) {
        for (const hit* firstHit : hitsByPlane[0]) {
            for (const hit* secondHit : hitsByPlane[1]) {
                vector<vector<const hit*>> candidates =
                    buildTracksFromSeedGenerous(
                        hitsByPlane, axis,
                        firstHit, secondHit,
                        cfg,
                        ambiguitySkipped);

                if (ambiguitySkipped > 0)
                    return {};

                foundTracks.insert(foundTracks.end(),
                                   candidates.begin(),
                                   candidates.end());
            }
        }
        return foundTracks;
    }

    // =====================================================
    // ===== GREEDY MODE ==================================
    // =====================================================
    while (true) {
        vector<const hit*> foundTrack;

        // Try all seeds
        for (const hit* firstHit : hitsByPlane[0]) {
            for (const hit* secondHit : hitsByPlane[1]) {
                vector<const hit*> candidate =
                    buildTrackFromSeed(hitsByPlane, axis, firstHit, secondHit, cfg, ambiguitySkipped);

                if (ambiguitySkipped > 0) {
                    return {};
                }

                if (!candidate.empty()) {
                    foundTrack = candidate;
                    break;
                }
            }
            if (!foundTrack.empty()) break;
        }

        if (foundTrack.empty()) {
            break;  // no more tracks
        }

        foundTracks.push_back(foundTrack);

        // Remove found hits (only if greedy is true)
        if (cfg.greedy) {
            for (size_t planeIndex = 0; planeIndex < foundTrack.size(); ++planeIndex) {
                vector<const hit*>& plane = hitsByPlane[planeIndex];
                plane.erase(remove(plane.begin(), plane.end(), foundTrack[planeIndex]), plane.end());
            }
        }
    }

    return foundTracks;
}

#endif // KALMAN_H