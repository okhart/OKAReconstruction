#include <TFile.h>
#include <TTree.h>
#include <TVector.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <TH2F.h>
#include "track.h"
#include "hough.h"
#include "filter.h"
#include "kalman.h"
#include "settings.h"

using namespace std;

// =====================================================
// ===== MAIN FUNCTION =================================
// =====================================================
int main() {

    // =====================================================
    // ===== OPEN INPUT FILE ===============================
    // =====================================================
    TFile *file = TFile::Open(INPUT_FILE, "READ");
    if (!file || file->IsZombie()) {
        cerr << "Error: Cannot open file " << INPUT_FILE << endl;
        return 1;
    }

    // =====================================================
    // ===== GET TREE AND BRANCHES =========================
    // =====================================================
    TTree *tree = (TTree*)file->Get(TREE_NAME);
    if (!tree) {
        cerr << "Error: Cannot find tree " << TREE_NAME << endl;
        file->Close();
        return 1;
    }

    TObjArray *branches = tree->GetListOfBranches();
    int nHits = branches->GetEntries() - 3;  // first 3 branches are truth arrays

    if (nHits < 2) {
        cerr << "Error: Need at least 2 hits" << endl;
        file->Close();
        return 1;
    }

    // =====================================================
    // ===== PREPARE BRANCH STORAGE =========================
    // =====================================================

    // Truth branches (only used in Simulation mode)
    float arr0[4], arr1[4], arr2[4];
    if (DATA_MODE == DataMode::Simulation) {
        ((TBranch*)branches->At(0))->SetAddress(arr0);
        ((TBranch*)branches->At(1))->SetAddress(arr1);
        ((TBranch*)branches->At(2))->SetAddress(arr2);
    }

    // Hit branches
    vector<vector<float>*> vecs(nHits, nullptr);
    for (int i = 0; i < nHits; i++) {
        TBranch* branch = (TBranch*)branches->At(i + 3);
        vecs[i] = nullptr;
        branch->SetAddress(&vecs[i]);
    }

    // =====================================================
    // ===== EVENT LOOP =====================================
    // =====================================================

    Long64_t nEntries = tree->GetEntries();
    if (MAX_EVENTS > 0 && MAX_EVENTS < nEntries) nEntries = MAX_EVENTS;

    cout << "Processing " << nEntries << " events..." << endl;

    int totalEvents = 0;
    int goodEvents = 0;
    int ambiguitySkippedEvents = 0;
    int lessTracksSkipped = 0;
    int chi2SkippedEvents = 0;

    for (Long64_t evt = 0; evt < nEntries; evt++) {

        // =====================================================
        // ===== READ EVENT ====================================
        // =====================================================
        tree->GetEntry(evt);

        // =====================================================
        // ===== BUILD HITS ====================================
        // =====================================================
        vector<hit*> trackHits;
        trackHits.reserve(nHits);

        for (int i = 0; i < nHits; i++) {
            vector<float>* vec = vecs[i];
            if (!vec || vec->size() < 4) continue;

            double z_val = (*vec)[0];
            double x_val = (*vec)[1];
            double y_val = (*vec)[2];
            double sigma_val = (*vec)[3];

            bool has_x = !std::isnan(x_val);
            bool has_y = !std::isnan(y_val);

            PlaneType type;
            if (has_x && !has_y) type = PlaneType::ZX;
            else if (!has_x && has_y) type = PlaneType::ZY;
            else type = PlaneType::ZXY;

            hit* h = new hit();
            h->setZ(z_val);
            h->setX(x_val);
            h->setY(y_val);
            h->setSigma(sigma_val);
            h->setType(type);
            trackHits.push_back(h);
        }

        // =====================================================
        // ===== CHECK ZXY AMBIGUITY ===========================
        // =====================================================
        // Skip event if there are identical ZXY hits

        bool ambiguousZXY = false;

        for (size_t i = 0; i < nHits; i++) {

            const hit* h1 = trackHits[i];

            if (h1->getType() != PlaneType::ZXY)
                continue;

            for (size_t j = i + 1; j < nHits; j++) {

                const hit* h2 = trackHits[j];

                if (h2->getType() != PlaneType::ZXY)
                    continue;

                bool sameZ =
                    h1->getZ() == h2->getZ();

                bool sameX =
                    h1->getX() == h2->getX();

                bool sameY =
                    h1->getY() == h2->getY();

                if (sameZ && sameX && sameY) {

                    ambiguousZXY = true;

                    break;
                }
                if (ambiguousZXY)
                    break;
            }

            if (ambiguousZXY)
                break;
        }

        if (ambiguousZXY){
            for (auto* h : trackHits) delete h;
            continue;
        }

        totalEvents++;

        // =====================================================
        // ===== RUN TRACK FINDER ==============================
        // =====================================================

        vector<vector<const hit*>> candidatesZX, candidatesZY;


        if (TRACKING_ALGORITHM == TrackingAlgorithm::Kalman) {
            // ---- Kalman filter ----
            CKFConfig ckfCfg;
            ckfCfg.minHits = CKF_MIN_HITS;
            ckfCfg.chi2Cut = CKF_CHI2_CUT;
            ckfCfg.processNoiseSlope = CKF_PROCESS_NOISE;
            ckfCfg.chi2AmbiguityCut = CKF_AMBIGUITY_CHI2;
            ckfCfg.setGreedy(CKF_GREEDY);
            ckfCfg.setGenerous(CKF_GENEROUS);

            auto hitsByPlaneZX = groupHitsByPlane(trackHits, PlaneType::ZX);
            auto hitsByPlaneZY = groupHitsByPlane(trackHits, PlaneType::ZY);

            int ambiguitySkippedThisEvent = 0;
            candidatesZX = findTracksCKF(hitsByPlaneZX, PlaneType::ZX,
                                         ambiguitySkippedThisEvent, ckfCfg);
            candidatesZY = findTracksCKF(hitsByPlaneZY, PlaneType::ZY,
                                         ambiguitySkippedThisEvent, ckfCfg);

            if (ambiguitySkippedThisEvent > 0) {
                for (auto* h : trackHits) delete h;
                candidatesZX.clear();
                candidatesZY.clear();
                ambiguitySkippedEvents++;
                continue;
            }
            
        } else {
            // ---- Hough transform ----
            HoughConfig hCfg;
            hCfg.a_min = HOUGH_A_MIN;
            hCfg.a_max = HOUGH_A_MAX;
            hCfg.n_a_bins = HOUGH_A_BINS;
            hCfg.n_b_bins = HOUGH_B_BINS;
            hCfg.threshold = HOUGH_THRESHOLD;
            hCfg.overlap = HOUGH_OVERLAP;
            hCfg.debug = DEBUG;        // use global DEBUG flag
            hCfg.b_min = HOUGH_B_MIN;
            hCfg.b_max = HOUGH_B_MAX;

            TFile* houghFile = nullptr;
            if (hCfg.debug) {
                houghFile = new TFile("hough_debug.root", "RECREATE");
            }

            vector<vector<HBin>> accZX, accZY;
            double bminZX, bmaxZX, bminZY, bmaxZY;

            buildHough(trackHits, PlaneType::ZX, hCfg, accZX, bminZX, bmaxZX,
                       houghFile, candidatesZX);
            buildHough(trackHits, PlaneType::ZY, hCfg, accZY, bminZY, bmaxZY,
                       houghFile, candidatesZY);

            if (houghFile) {
                houghFile->Close();
                delete houghFile;
            }
        }

        // Check if we have enough tracks in both projections
        if ((int)candidatesZX.size() < EXPECTED_N_TRACKS ||
            (int)candidatesZY.size() < EXPECTED_N_TRACKS) {
            lessTracksSkipped++;
            for (auto* h : trackHits) delete h;
            candidatesZX.clear();
            candidatesZY.clear();
            continue;
        }

        // =====================================================
        // ===== COMBINE ZX AND ZY CANDIDATES ==================
        // =====================================================

        vector<track*> allTracks;

        for (auto& zx_hits : candidatesZX) {
            for (auto& zy_hits : candidatesZY) {
                // Find common ZXY hit
                const hit* common = nullptr;
                for (auto* h1 : zx_hits) {
                    if (h1->getType() != PlaneType::ZXY) continue;
                    for (auto* h2 : zy_hits) {
                        if (h2->getType() != PlaneType::ZXY) continue;
                        if (h1 == h2) {
                            common = h1;
                            break;
                        }
                    }
                    if (common) break;
                }
                if (!common) continue;

                // Build combined hit list
                vector<const hit*> all;
                all.push_back(common);
                for (auto* h : zx_hits) if (h != common) all.push_back(h);
                for (auto* h : zy_hits) if (h != common) all.push_back(h);

                const hit** hitArray = new const hit*[all.size()];
                for (size_t i = 0; i < all.size(); i++) hitArray[i] = all[i];
                track* t = new track(hitArray, all.size());
                delete[] hitArray;
                allTracks.push_back(t);
            }
        }

        // =====================================================
        // ===== EXPAND / FILTER TRACKS =======================
        // =====================================================

        vector<track*> validTracks;

        // If using Kalman, each candidate already has exactly one hit per plane,
        // so we just take them as they are.
        if (TRACKING_ALGORITHM == TrackingAlgorithm::Kalman) {
            for (auto* cand : allTracks) {
                if (cand->getNHits() != 5) continue;  // should always be 5
                const hit** hitArray = new const hit*[5];
                for (int i = 0; i < 5; i++) hitArray[i] = cand->getHit(i);
                track* t = new track(hitArray, 5);
                delete[] hitArray;
                validTracks.push_back(t);
            }
            for (auto* t : allTracks) delete t;
            allTracks.clear();
        } else {
            // For Hough, we need to build all 5-hit combinations (2 ZX + 2 ZY + 1 ZXY)
            for (auto* cand : allTracks) {
                vector<const hit*> zxHits, zyHits, zxyHits;
                for (int i = 0; i < cand->getNHits(); i++) {
                    const hit* h = cand->getHit(i);
                    if (h->getType() == PlaneType::ZX) zxHits.push_back(h);
                    else if (h->getType() == PlaneType::ZY) zyHits.push_back(h);
                    else if (h->getType() == PlaneType::ZXY) zxyHits.push_back(h);
                }
                if (zxHits.size() < 2 || zyHits.size() < 2 || zxyHits.size() < 1) continue;

                for (size_t iz = 0; iz < zxyHits.size(); iz++) {
                    for (size_t i = 0; i < zxHits.size(); i++) {
                        for (size_t j = i+1; j < zxHits.size(); j++) {
                            for (size_t k = 0; k < zyHits.size(); k++) {
                                for (size_t l = k+1; l < zyHits.size(); l++) {
                                    const hit* arr[5] = {
                                        zxyHits[iz],
                                        zxHits[i], zxHits[j],
                                        zyHits[k], zyHits[l]
                                    };
                                    track* t = new track(arr, 5);
                                    validTracks.push_back(t);
                                }
                            }
                        }
                    }
                }
            }
            for (auto* t : allTracks) delete t;
            allTracks.clear();
        }
        
        vector<track*> finalTracks;

        // For generous Kalman, keep only the best EXPECTED_N_TRACKS tracks
        if (TRACKING_ALGORITHM == TrackingAlgorithm::Kalman && CKF_GENEROUS) {
            finalTracks = filterTracks(validTracks);
            if ((int)finalTracks.size() > EXPECTED_N_TRACKS) {
                finalTracks.resize(EXPECTED_N_TRACKS);
            }
        }
        else if (TRACKING_ALGORITHM == TrackingAlgorithm::Hough){
            finalTracks = filterTracks(validTracks);
        }
        else{
            finalTracks = validTracks;
        }

        // =====================================================
        // ===== APPLY QUALITY CUTS ============================
        // =====================================================

        bool passedQualityCuts = true;

        // Apply χ² cut to all tracks
        for (auto* reco : finalTracks) {
            if (reco->getChi2() > CHI2_PER_NDF_CUT) {
                chi2SkippedEvents++;
                passedQualityCuts = false;
                break;
            }
        }

        if (!passedQualityCuts) {
            for (auto* t : finalTracks) delete t;
            for (auto* h : trackHits) delete h;
            candidatesZX.clear();
            candidatesZY.clear();
            continue;
        }

        // =====================================================
        // ===== EVALUATE (depends on DataMode) ================
        // =====================================================

        if (DATA_MODE == DataMode::Simulation) {
            // Compare with truth
            struct TrueTrack { double a,b,c,d; };
            vector<TrueTrack> truth = {
                {arr0[0], arr0[1], arr0[2], arr0[3]},
                {arr1[0], arr1[1], arr1[2], arr1[3]},
                {arr2[0], arr2[1], arr2[2], arr2[3]}
            };

            vector<bool> truthUsed(EXPECTED_N_TRACKS, false);
            int matchedTracks = 0;

            for (auto* reco : finalTracks) {
                for (int t = 0; t < EXPECTED_N_TRACKS; t++) {
                    if (truthUsed[t]) continue;
                    int match = 0;
                    // Each true track has 5 hits (indices 0..4)
                    for (int ib = 0; ib < 5; ib++) {
                        hit* th = trackHits[t*5 + ib];
                        for (int i = 0; i < reco->getNHits(); i++) {
                            if (reco->getHit(i) == th) { match++; break; }
                        }
                    }
                    if (match == 5) {
                        truthUsed[t] = true;
                        matchedTracks++;
                        break;
                    }
                }
            }

            if (matchedTracks == EXPECTED_N_TRACKS) {
                goodEvents++;
            } else if (DEBUG) {
                // --- Debug output for incomplete matching ---
                cout << "\n=== Incomplete match: matchedTracks = "
                     << matchedTracks << "/" << EXPECTED_N_TRACKS << " ===" << endl;

                // Print true tracks and their hits
                cout << "--- True tracks ---" << endl;
                for (int t = 0; t < EXPECTED_N_TRACKS; t++) {
                    cout << "True track " << t
                         << ": a=" << truth[t].a
                         << ", b=" << truth[t].b
                         << ", c=" << truth[t].c
                         << ", d=" << truth[t].d << endl;
                    for (int k = 0; k < 5; k++) {
                        hit* th = trackHits[t*5 + k];
                        cout << "    Hit: Z=" << th->getZ();
                        if (th->hasX()) cout << ", X=" << th->getX();
                        if (th->hasY()) cout << ", Y=" << th->getY();
                        cout << ", Sigma=" << th->getSigma() << endl;
                    }
                }

                // Print reconstructed tracks
                cout << "--- Reconstructed tracks ---" << endl;
                for (size_t r = 0; r < finalTracks.size(); r++) {
                    cout << "Reco track " << r << ":" << endl;
                    finalTracks[r]->print();
                }
                cout << "=======================================" << endl;
            }
        } else {
            // Experimental: just count events that passed all cuts
            goodEvents++;
        }

        // =====================================================
        // ===== MEMORY CLEANUP ================================
        // =====================================================
        for (auto* t : finalTracks) delete t;
        candidatesZX.clear();
        candidatesZY.clear();
        for (auto* h : trackHits) delete h;
    }

    // =====================================================
    // ===== FINAL STATISTICS ==============================
    // =====================================================

    if (PRINT_STATISTICS) {
        cout << "\n==============================" << endl;
        cout << "TOTAL EVENTS PROCESSED: " << totalEvents << endl;
        cout << "SKIPPED (ambiguous hits): " << ambiguitySkippedEvents << endl;
        cout << "SKIPPED (< " << EXPECTED_N_TRACKS << " tracks): " << lessTracksSkipped << endl;
        cout << "SKIPPED (χ² cut):         " << chi2SkippedEvents << endl;

        if (DATA_MODE == DataMode::Simulation) {
            cout << "GOOD EVENTS (fully matched): " << goodEvents << endl;
            cout << "EFFICIENCY:                  "
                 << (double)(goodEvents + ambiguitySkippedEvents + chi2SkippedEvents + lessTracksSkipped)/ totalEvents * 100 << "%" << endl;
        } else {
            cout << "GOOD EVENTS (passed all cuts): " << goodEvents << endl;
            cout << "ACCEPTANCE:                    "
                 << (double)goodEvents / totalEvents * 100 << "%" << endl;
        }
        cout << "==============================" << endl;
    }

    file->Close();
    cout << "\nDone." << endl;
    return 0;
}
