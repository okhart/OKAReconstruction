#ifndef HOUGH_H
#define HOUGH_H

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

static int g_hough_counter = 0;  // Global counter for debug output

// =====================================================
// ===== HOUGH BIN STRUCTURE ===========================
// =====================================================

/**
 * Represents a single bin in the Hough accumulator.
 * Stores the number of hits and pointers to those hits.
 */
struct HBin {
    int count = 0;
    vector<const hit*> hits;
};

// =====================================================
// ===== HOUGH CONFIGURATION STRUCTURE =================
// =====================================================

/**
 * Configuration parameters for the Hough transform.
 * These are normally filled from the global settings.
 */
struct HoughConfig {
    double a_min;       ///< Minimum slope
    double a_max;       ///< Maximum slope
    int n_a_bins;       ///< Number of bins in slope dimension
    int n_b_bins;       ///< Number of bins in intercept dimension
    int threshold;      ///< Minimum hits to consider a peak
    double overlap;     ///< Fractional overlap between bins
    bool debug;         ///< Enable debug output and histograms
    double b_min;       ///< Manually fixed intercept min (if both zero, auto-calc)
    double b_max;       ///< Manually fixed intercept max
};

// =====================================================
// ===== BUILD HOUGH SPACE =============================
// =====================================================

/**
 * Performs the Hough transform to find candidate lines in one projection (ZX or ZY).
 *
 * @param hits          All hits from the event.
 * @param axis          Which projection (PlaneType::ZX or PlaneType::ZY).
 * @param cfg           Hough configuration.
 * @param acc           Output accumulator grid (2D vector of HBin).
 * @param b_min         Output: minimum intercept found (or fixed).
 * @param b_max         Output: maximum intercept found (or fixed).
 * @param outFile       ROOT file for debug histograms (if debug enabled).
 * @param outCandidates Output: vector of hit groups, each group corresponds to a line candidate.
 */
void buildHough(
    const vector<hit*>& hits,
    PlaneType axis,
    const HoughConfig& cfg,
    vector<vector<HBin>>& acc,
    double& b_min,
    double& b_max,
    TFile* outFile,
    vector<vector<const hit*>>& outCandidates)
{
    const double a_min = cfg.a_min;
    const double a_max = cfg.a_max;
    const double da = (a_max - a_min) / cfg.n_a_bins;

    // =====================================================
    // ===== DETERMINE B RANGE ============================
    // =====================================================

    if (cfg.b_min != 0.0 || cfg.b_max != 0.0) {
        // Use fixed values provided by the user
        b_min = cfg.b_min;
        b_max = cfg.b_max;
        if (cfg.debug) {
            cout << "[Hough] Using fixed b range: [" << b_min << ", " << b_max << "]" << endl;
        }
    } else {
        // Auto-calculate from data
        b_min = std::numeric_limits<double>::max();
        b_max = std::numeric_limits<double>::lowest();

        for (const hit* h : hits) {
            if (!h->hasCoord(axis)) continue;
            double coord = h->getCoord(axis);
            double z = h->getZ();

            for (int ia = 0; ia < cfg.n_a_bins; ia++) {
                double a = cfg.a_min + ia * da;
                double b = coord - a * z;
                b_min = min(b_min, b);
                b_max = max(b_max, b);
            }
        }

        // Add a small margin
        double margin = 0.1 * (b_max - b_min);
        b_min -= margin;
        b_max += margin;

        if (cfg.debug) {
            cout << "[Hough] Auto-calculated b range: [" << b_min << ", " << b_max << "]" << endl;
        }
    }

    const double db = (b_max - b_min) / cfg.n_b_bins;

    // Initialize accumulator
    acc.assign(cfg.n_a_bins, vector<HBin>(cfg.n_b_bins));

    // =====================================================
    // ===== BIN-FIRST LOOP ================================
    // =====================================================

    for (int ia = 0; ia < cfg.n_a_bins; ia++) {
        double a0 = a_min + ia * da;
        double a1 = a0 + da;

        for (int ib = 0; ib < cfg.n_b_bins; ib++) {
            double b0 = b_min + ib * db;
            double b1 = b0 + db;

            // Apply overlap margin
            double a_margin = cfg.overlap * da;
            double b_margin = cfg.overlap * db;
            b0 -= b_margin;
            b1 += b_margin;

            HBin& bin = acc[ia][ib];

            // Check each hit
            for (const hit* h : hits) {
                if (!h->hasCoord(axis)) continue;

                double z = h->getZ();
                double coord = h->getCoord(axis);

                double b_at_a0 = coord - a0 * z;
                double b_at_a1 = coord - a1 * z;
                double b_min_line = min(b_at_a0, b_at_a1);
                double b_max_line = max(b_at_a0, b_at_a1);

                if (b_max_line < b0 || b_min_line > b1)
                    continue;

                bin.hits.push_back(h);
                bin.count++;
            }

            // ===== PEAK DETECTION =========================
            if (bin.count >= cfg.threshold) {
                if (cfg.debug) {
                    cout << ++g_hough_counter << "-" << ia << " " << ib << " " << bin.count << endl;
                }

                vector<const hit*> ordered;
                ordered.reserve(bin.count);

                // Prioritise ZXY hits (they carry both coordinates)
                for (auto* h : bin.hits)
                    if (h->getType() == PlaneType::ZXY)
                        ordered.push_back(h);

                // Then single-coordinate hits
                for (auto* h : bin.hits)
                    if (h->getType() != PlaneType::ZXY)
                        ordered.push_back(h);

                outCandidates.push_back(ordered);
            }
        }
    }

    // =====================================================
    // ===== DEBUG HISTOGRAM ===============================
    // =====================================================

    if (cfg.debug) {
        static int hough_id = 0;
        TString name = Form("hough_%d", hough_id++);
        TH2F* h2 = new TH2F(
            name,
            "Hough space; a; b",
            cfg.n_a_bins, cfg.a_min, cfg.a_max,
            cfg.n_b_bins, b_min, b_max
        );

        for (int ia = 0; ia < cfg.n_a_bins; ia++) {
            for (int ib = 0; ib < cfg.n_b_bins; ib++) {
                if (acc[ia][ib].count > 0) {
                    h2->SetBinContent(ia + 1, ib + 1, acc[ia][ib].count);
                }
            }
        }

        outFile->cd();
        h2->Write();
        delete h2;
        cout << "[Hough] Histogram saved: " << name << endl;
    }
}

#endif // HOUGH_H