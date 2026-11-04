#include <TFile.h>
#include <TTree.h>
#include <TVector.h>

#include <iostream>
#include <cmath>
#include <vector>

#include "track.h"
#include "settings.h"
#include "tracking.h"

using namespace std;

int main() {
    // -------------------------------------------------
    // Open input ROOT file and get the event tree.
    // -------------------------------------------------
    TFile* file = TFile::Open(INPUT_FILE, "READ");
    if (!file || file->IsZombie()) {
        cerr << "Error: Cannot open file " << INPUT_FILE << endl;
        return 1;
    }

    TTree* tree = static_cast<TTree*>(file->Get(TREE_NAME));
    if (!tree) {
        cerr << "Error: Cannot find tree " << TREE_NAME << endl;
        file->Close();
        return 1;
    }

    TObjArray* branches = tree->GetListOfBranches();
    const int nHitBranches = branches->GetEntries() - 3;
    if (nHitBranches < 2) {
        cerr << "Error: Need at least 2 hit branches" << endl;
        file->Close();
        return 1;
    }

    // -------------------------------------------------
    // Set addresses of truth and hit branches.
    // -------------------------------------------------
    float arr0[4] = {0}, arr1[4] = {0}, arr2[4] = {0};
    if (DATA_MODE == DataMode::Simulation) {
        static_cast<TBranch*>(branches->At(0))->SetAddress(arr0);
        static_cast<TBranch*>(branches->At(1))->SetAddress(arr1);
        static_cast<TBranch*>(branches->At(2))->SetAddress(arr2);
    }

    vector<vector<float>*> hitBranchData(nHitBranches, nullptr);
    for (int i = 0; i < nHitBranches; ++i) {
        TBranch* branch = static_cast<TBranch*>(branches->At(i + 3));
        branch->SetAddress(&hitBranchData[i]);
    }

    // -------------------------------------------------
    // Event loop.
    // -------------------------------------------------
    Long64_t nEntries = tree->GetEntries();
    if (MAX_EVENTS > 0 && MAX_EVENTS < nEntries)
        nEntries = MAX_EVENTS;

    cout << "Processing " << nEntries << " events..." << endl;

    int totalEvents = 0;
    int goodEvents = 0;
    int ambiguitySkippedEvents = 0;
    int lessTracksSkipped = 0;
    int chi2SkippedEvents = 0;

    for (Long64_t evt = 0; evt < nEntries; ++evt) {
        tree->GetEntry(evt);

        // -------------------------------------------------
        // Build hit objects for this event.
        // -------------------------------------------------
        vector<hit*> eventHits;
        eventHits.reserve(nHitBranches);

        for (int i = 0; i < nHitBranches; ++i) {
            vector<float>* values = hitBranchData[i];
            if (!values || values->size() < 4)
                continue;

            const double z = (*values)[0];
            const double x = (*values)[1];
            const double y = (*values)[2];
            const double sigma = (*values)[3];

            const bool hasX = !std::isnan(x);
            const bool hasY = !std::isnan(y);

            PlaneType type;
            if (hasX && !hasY)
                type = PlaneType::ZX;
            else if (!hasX && hasY)
                type = PlaneType::ZY;
            else
                type = PlaneType::ZXY;

            hit* h = new hit();
            h->setZ(z);
            h->setX(x);
            h->setY(y);
            h->setSigma(sigma);
            h->setType(type);
            eventHits.push_back(h);
        }

        // -------------------------------------------------
        // Reject events with duplicated ZXY hits.
        // -------------------------------------------------
        bool ambiguousZXY = false;
        for (size_t i = 0; i < eventHits.size() && !ambiguousZXY; ++i) {
            const hit* h1 = eventHits[i];
            if (h1->getType() != PlaneType::ZXY)
                continue;

            for (size_t j = i + 1; j < eventHits.size(); ++j) {
                const hit* h2 = eventHits[j];
                if (h2->getType() != PlaneType::ZXY)
                    continue;

                if (h1->getZ() == h2->getZ() &&
                    h1->getX() == h2->getX() &&
                    h1->getY() == h2->getY()) {
                    ambiguousZXY = true;
                    break;
                }
            }
        }

        if (ambiguousZXY) {
            for (hit* h : eventHits)
                delete h;
            continue;
        }

        ++totalEvents;

        // -------------------------------------------------
        // Reconstruct straight tracks in the supplied hit set.
        // This function is independent of the magnet. It can later be
        // called separately for hitsBefore and hitsAfter.
        // -------------------------------------------------
        int ambiguitySkippedThisEvent = 0;
        vector<track> finalTracks = findTracks(eventHits, ambiguitySkippedThisEvent);

        if (ambiguitySkippedThisEvent > 0) {
            ++ambiguitySkippedEvents;
            for (hit* h : eventHits)
                delete h;
            continue;
        }

        if (static_cast<int>(finalTracks.size()) < EXPECTED_N_TRACKS) {
            ++lessTracksSkipped;
            for (hit* h : eventHits)
                delete h;
            continue;
        }

        // Keep only the best expected number in generous CKF mode.
        if (TRACKING_ALGORITHM == TrackingAlgorithm::Kalman && CKF_GENEROUS &&
            static_cast<int>(finalTracks.size()) > EXPECTED_N_TRACKS) {
            finalTracks.resize(EXPECTED_N_TRACKS);
        }

        // -------------------------------------------------
        // Apply track quality cuts.
        // -------------------------------------------------
        bool passedQualityCuts = true;
        for (const track& reconstructed : finalTracks) {
            if (reconstructed.getChi2() > CHI2_PER_NDF_CUT) {
                ++chi2SkippedEvents;
                passedQualityCuts = false;
                break;
            }
        }

        if (!passedQualityCuts) {
            for (hit* h : eventHits)
                delete h;
            continue;
        }

        // -------------------------------------------------
        // Evaluate reconstruction (simulation or experimental mode).
        // -------------------------------------------------
        if (DATA_MODE == DataMode::Simulation) {
            struct TrueTrack {
                double a, b, c, d;
            };

            vector<TrueTrack> truth = {
                {arr0[0], arr0[1], arr0[2], arr0[3]},
                {arr1[0], arr1[1], arr1[2], arr1[3]},
                {arr2[0], arr2[1], arr2[2], arr2[3]}
            };

            vector<bool> truthUsed(EXPECTED_N_TRACKS, false);
            int matchedTracks = 0;

            for (const track& reconstructed : finalTracks) {
                for (int truthIndex = 0; truthIndex < EXPECTED_N_TRACKS; ++truthIndex) {
                    if (truthUsed[truthIndex])
                        continue;

                    int matchedHits = 0;
                    for (int hitIndex = 0; hitIndex < 5; ++hitIndex) {
                        const size_t eventHitIndex =
                            static_cast<size_t>(truthIndex * 5 + hitIndex);
                        if (eventHitIndex >= eventHits.size())
                            continue;

                        const hit* truthHit = eventHits[eventHitIndex];
                        for (int i = 0; i < reconstructed.getNHits(); ++i) {
                            if (reconstructed.getHit(i) == truthHit) {
                                ++matchedHits;
                                break;
                            }
                        }
                    }

                    if (matchedHits == 5) {
                        truthUsed[truthIndex] = true;
                        ++matchedTracks;
                        break;
                    }
                }
            }

            if (matchedTracks == EXPECTED_N_TRACKS) {
                ++goodEvents;
            } else if (DEBUG) {
                cout << "\n=== Incomplete match: matchedTracks = "
                     << matchedTracks << "/" << EXPECTED_N_TRACKS << " ===" << endl;

                cout << "--- True tracks ---" << endl;
                for (int t = 0; t < EXPECTED_N_TRACKS; ++t) {
                    cout << "True track " << t
                         << ": a=" << truth[t].a
                         << ", b=" << truth[t].b
                         << ", c=" << truth[t].c
                         << ", d=" << truth[t].d << endl;

                    for (int k = 0; k < 5; ++k) {
                        const size_t index = static_cast<size_t>(t * 5 + k);
                        if (index >= eventHits.size())
                            continue;
                        const hit* h = eventHits[index];
                        cout << "    Hit: Z=" << h->getZ();
                        if (h->hasX()) cout << ", X=" << h->getX();
                        if (h->hasY()) cout << ", Y=" << h->getY();
                        cout << ", Sigma=" << h->getSigma() << endl;
                    }
                }

                cout << "--- Reconstructed tracks ---" << endl;
                for (size_t i = 0; i < finalTracks.size(); ++i) {
                    cout << "Reco track " << i << ":" << endl;
                    finalTracks[i].print();
                }
                cout << "=======================================" << endl;
            }
        } else {
            ++goodEvents;
        }

        // track objects own their vector of hit pointers by value;
        // the hit objects themselves are owned and deleted here.
        for (hit* h : eventHits)
            delete h;
    }

    // -------------------------------------------------
    // Final statistics.
    // -------------------------------------------------
    if (PRINT_STATISTICS) {
        cout << "\n==============================" << endl;
        cout << "TOTAL EVENTS PROCESSED: " << totalEvents << endl;
        cout << "SKIPPED (ambiguous hits): " << ambiguitySkippedEvents << endl;
        cout << "SKIPPED (< " << EXPECTED_N_TRACKS << " tracks): "
             << lessTracksSkipped << endl;
        cout << "SKIPPED (chi2 cut): " << chi2SkippedEvents << endl;

        if (DATA_MODE == DataMode::Simulation) {
            cout << "GOOD EVENTS (fully matched): " << goodEvents << endl;
            if (totalEvents > 0) {
                cout << "EFFICIENCY: "
                     << static_cast<double>(goodEvents + ambiguitySkippedEvents +
                                            chi2SkippedEvents + lessTracksSkipped)
                        / totalEvents * 100.0
                     << "%" << endl;
            }
        } else {
            cout << "GOOD EVENTS (passed all cuts): " << goodEvents << endl;
            if (totalEvents > 0) {
                cout << "ACCEPTANCE: "
                     << static_cast<double>(goodEvents) / totalEvents * 100.0
                     << "%" << endl;
            }
        }
        cout << "==============================" << endl;
    }

    file->Close();
    cout << "\nDone." << endl;
    return 0;
}
