#ifndef TRACKING_H
#define TRACKING_H

#include <vector>

#include "track.h"
#include "settings.h"
#include "kalman.h"
#include "hough.h"
#include "filter.h"

using namespace std;

/**
 * Reconstructs straight tracks in one detector region.
 *
 * Input:
 *   hits - hits belonging to this region only.
 *
 * Output:
 *   reconstructed straight tracks. Each track contains its hit pointers,
 *   fitted ZX/ZY parameters and chi2.
 *
 * This function does not know about the magnet and does not bend tracks.
 * It can be called independently for the detector regions before and after
 * the magnet.
 *
 * ambiguitySkipped counts CKF candidates rejected due to ambiguity.
 */
inline vector<track> findTracks(
    const vector<Hit*>& hits,
    int& ambiguitySkipped)
{
    vector<vector<const Hit*> > candidatesZX;
    vector<vector<const Hit*> > candidatesZY;

    // -------------------------------------------------
    // Find projection candidates.
    // -------------------------------------------------
    if (TRACKING_ALGORITHM == TrackingAlgorithm::Kalman) {
        CKFConfig ckfCfg;
        ckfCfg.minHits = CKF_MIN_HITS;
        ckfCfg.chi2Cut = CKF_CHI2_CUT;
        ckfCfg.processNoiseSlope = CKF_PROCESS_NOISE;
        ckfCfg.chi2AmbiguityCut = CKF_AMBIGUITY_CHI2;
        ckfCfg.setGreedy(CKF_GREEDY);
        ckfCfg.setGenerous(CKF_GENEROUS);
        ckfCfg.zOrder = CKF_Z_ORDER;

        const vector<vector<const Hit*> > hitsByPlaneZX =
            groupHitsByPlane(hits, PlaneType::ZX, ckfCfg.zOrder);
        const vector<vector<const Hit*> > hitsByPlaneZY =
            groupHitsByPlane(hits, PlaneType::ZY, ckfCfg.zOrder);

        candidatesZX = findTracksCKF(
            hitsByPlaneZX, PlaneType::ZX, ambiguitySkipped, ckfCfg);
        candidatesZY = findTracksCKF(
            hitsByPlaneZY, PlaneType::ZY, ambiguitySkipped, ckfCfg);

        if (ambiguitySkipped > 0)
            return {};
    } else {
        HoughConfig hCfg;
        hCfg.a_min = HOUGH_A_MIN;
        hCfg.a_max = HOUGH_A_MAX;
        hCfg.n_a_bins = HOUGH_A_BINS;
        hCfg.n_b_bins = HOUGH_B_BINS;
        hCfg.threshold = HOUGH_THRESHOLD;
        hCfg.overlap = HOUGH_OVERLAP;
        hCfg.debug = DEBUG;
        hCfg.b_min = HOUGH_B_MIN;
        hCfg.b_max = HOUGH_B_MAX;

        TFile* houghFile = nullptr;
        if (hCfg.debug)
            houghFile = new TFile("hough_debug.root", "RECREATE");

        vector<vector<HBin> > accumulatorZX;
        vector<vector<HBin> > accumulatorZY;
        double bMinZX = 0.0, bMaxZX = 0.0;
        double bMinZY = 0.0, bMaxZY = 0.0;

        buildHough(hits, PlaneType::ZX, hCfg,
                   accumulatorZX, bMinZX, bMaxZX,
                   houghFile, candidatesZX);
        buildHough(hits, PlaneType::ZY, hCfg,
                   accumulatorZY, bMinZY, bMaxZY,
                   houghFile, candidatesZY);

        if (houghFile) {
            houghFile->Close();
            delete houghFile;
        }
    }

    // The original event-level reconstruction requires at least the
    // expected number of projection candidates in both views.
    if (static_cast<int>(candidatesZX.size()) < EXPECTED_N_TRACKS ||
        static_cast<int>(candidatesZY.size()) < EXPECTED_N_TRACKS)
        return {};

    // -------------------------------------------------
    // Match ZX and ZY candidates by their common ZXY hit.
    // -------------------------------------------------
    vector<track> combinedTracks;

    for (const vector<const Hit*>& zxHits : candidatesZX) {
        for (const vector<const Hit*>& zyHits : candidatesZY) {
            const Hit* commonHit = nullptr;

            for (const Hit* zxHit : zxHits) {
                if (zxHit->getType() != PlaneType::ZXY)
                    continue;

                for (const Hit* zyHit : zyHits) {
                    if (zyHit->getType() != PlaneType::ZXY)
                        continue;

                    if (zxHit == zyHit) {
                        commonHit = zxHit;
                        break;
                    }
                }

                if (commonHit)
                    break;
            }

            if (!commonHit)
                continue;

            vector<const Hit*> combinedHits;
            combinedHits.reserve(zxHits.size() + zyHits.size() - 1);
            combinedHits.push_back(commonHit);

            for (const Hit* h : zxHits) {
                if (h != commonHit)
                    combinedHits.push_back(h);
            }

            for (const Hit* h : zyHits) {
                if (h != commonHit)
                    combinedHits.push_back(h);
            }

            combinedTracks.push_back(track(combinedHits));
        }
    }

    // -------------------------------------------------
    // Build the final 3D track candidates.
    // -------------------------------------------------
    vector<track> validTracks;

    if (TRACKING_ALGORITHM == TrackingAlgorithm::Kalman) {
        // CKF candidates already contain one hit per plane.
        for (const track& candidate : combinedTracks) {
            if (candidate.getNHits() == 5)
                validTracks.push_back(candidate);
        }
    } else {
        // Hough candidates may contain several hits in each projection.
        // Build all combinations with 2 ZX, 2 ZY and 1 common ZXY hit.
        for (const track& candidate : combinedTracks) {
            vector<const Hit*> zxHits;
            vector<const Hit*> zyHits;
            vector<const Hit*> zxyHits;

            for (int i = 0; i < candidate.getNHits(); ++i) {
                const Hit* h = candidate.getHit(i);
                if (h->getType() == PlaneType::ZX)
                    zxHits.push_back(h);
                else if (h->getType() == PlaneType::ZY)
                    zyHits.push_back(h);
                else if (h->getType() == PlaneType::ZXY)
                    zxyHits.push_back(h);
            }

            for (const Hit* zxyHit : zxyHits) {
                for (size_t i = 0; i < zxHits.size(); ++i) {
                    for (size_t j = i + 1; j < zxHits.size(); ++j) {
                        for (size_t k = 0; k < zyHits.size(); ++k) {
                            for (size_t l = k + 1; l < zyHits.size(); ++l) {
                                vector<const Hit*> selectedHits = {
                                    zxyHit,
                                    zxHits[i], zxHits[j],
                                    zyHits[k], zyHits[l]
                                };
                                validTracks.push_back(track(selectedHits));
                            }
                        }
                    }
                }
            }
        }
    }

    // -------------------------------------------------
    // Remove overlaps in the same modes as before.
    // -------------------------------------------------
    if ((TRACKING_ALGORITHM == TrackingAlgorithm::Kalman && CKF_GENEROUS) ||
        TRACKING_ALGORITHM == TrackingAlgorithm::Hough) {
        return filterTracks(validTracks);
    }

    return validTracks;
}

// Simple interface when the caller does not need the ambiguity counter.
inline vector<track> findTracks(const vector<Hit*>& hits)
{
    int ambiguitySkipped = 0;
    return findTracks(hits, ambiguitySkipped);
}

#endif // TRACKING_H
