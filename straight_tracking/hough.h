#ifndef HOUGH_H
#define HOUGH_H

#include <TFile.h>
#include <TTree.h>
#include <TVector.h>
#include <TH2F.h>

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>

#include "track.h"

using namespace std;

static int g_hough_counter = 0;

// =====================================================
// ===== HOUGH BIN STRUCTURE ===========================
// =====================================================
struct HBin {
    int count = 0;
    vector<const hit*> hits;
};

// =====================================================
// ===== HOUGH CONFIGURATION STRUCTURE =================
// =====================================================
struct HoughConfig {
    double a_min;
    double a_max;
    int n_a_bins;
    int n_b_bins;
    int threshold;
    double overlap;
    bool debug;
    double b_min;
    double b_max;
};

// =====================================================
// ===== BUILD HOUGH SPACE =============================
// =====================================================
/**
 * Performs the Hough transform in one projection (ZX or ZY).
 */
inline void buildHough(
    const vector<hit*>& hits,
    PlaneType axis,
    const HoughConfig& cfg,
    vector<vector<HBin> >& acc,
    double& b_min,
    double& b_max,
    TFile* outFile,
    vector<vector<const hit*> >& outCandidates)
{
    const double a_min = cfg.a_min;
    const double a_max = cfg.a_max;
    const double da = (a_max - a_min) / cfg.n_a_bins;

    // Determine intercept range.
    if (cfg.b_min != 0.0 || cfg.b_max != 0.0) {
        b_min = cfg.b_min;
        b_max = cfg.b_max;
        if (cfg.debug) {
            cout << "[Hough] Using fixed b range: ["
                 << b_min << ", " << b_max << "]" << endl;
        }
    } else {
        b_min = numeric_limits<double>::max();
        b_max = numeric_limits<double>::lowest();

        for (const hit* h : hits) {
            if (!h->hasCoord(axis))
                continue;

            const double coord = h->getCoord(axis);
            const double z = h->getZ();

            for (int ia = 0; ia < cfg.n_a_bins; ++ia) {
                const double a = cfg.a_min + ia * da;
                const double b = coord - a * z;
                b_min = min(b_min, b);
                b_max = max(b_max, b);
            }
        }

        const double margin = 0.1 * (b_max - b_min);
        b_min -= margin;
        b_max += margin;

        if (cfg.debug) {
            cout << "[Hough] Auto-calculated b range: ["
                 << b_min << ", " << b_max << "]" << endl;
        }
    }

    const double db = (b_max - b_min) / cfg.n_b_bins;
    acc.assign(cfg.n_a_bins, vector<HBin>(cfg.n_b_bins));

    // Fill accumulator.
    for (int ia = 0; ia < cfg.n_a_bins; ++ia) {
        const double a0 = a_min + ia * da;
        const double a1 = a0 + da;

        for (int ib = 0; ib < cfg.n_b_bins; ++ib) {
            double b0 = b_min + ib * db;
            double b1 = b0 + db;
            const double a_margin = cfg.overlap * da;
            const double b_margin = cfg.overlap * db;
            b0 -= b_margin;
            b1 += b_margin;

            (void)a_margin;
            HBin& bin = acc[ia][ib];

            for (const hit* h : hits) {
                if (!h->hasCoord(axis))
                    continue;

                const double z = h->getZ();
                const double coord = h->getCoord(axis);
                const double b_at_a0 = coord - a0 * z;
                const double b_at_a1 = coord - a1 * z;
                const double b_min_line = min(b_at_a0, b_at_a1);
                const double b_max_line = max(b_at_a0, b_at_a1);

                if (b_max_line < b0 || b_min_line > b1)
                    continue;

                bin.hits.push_back(h);
                ++bin.count;
            }

            if (bin.count >= cfg.threshold) {
                if (cfg.debug)
                    cout << ++g_hough_counter << "-" << ia << " "
                         << ib << " " << bin.count << endl;

                vector<const hit*> ordered;
                ordered.reserve(bin.count);

                // Prioritise ZXY hits (they carry both coordinates).
                for (const hit* h : bin.hits) {
                    if (h->getType() == PlaneType::ZXY)
                        ordered.push_back(h);
                }

                // Then single-coordinate hits.
                for (const hit* h : bin.hits) {
                    if (h->getType() != PlaneType::ZXY)
                        ordered.push_back(h);
                }

                outCandidates.push_back(ordered);
            }
        }
    }

    // Save debug histogram.
    if (cfg.debug) {
        static int hough_id = 0;
        TString name = Form("hough_%d", hough_id++);
        TH2F* h2 = new TH2F(
            name,
            "Hough space; a; b",
            cfg.n_a_bins, cfg.a_min, cfg.a_max,
            cfg.n_b_bins, b_min, b_max);

        for (int ia = 0; ia < cfg.n_a_bins; ++ia) {
            for (int ib = 0; ib < cfg.n_b_bins; ++ib) {
                if (acc[ia][ib].count > 0)
                    h2->SetBinContent(ia + 1, ib + 1, acc[ia][ib].count);
            }
        }

        outFile->cd();
        h2->Write();
        delete h2;
        cout << "[Hough] Histogram saved: " << name << endl;
    }
}

#endif // HOUGH_H
