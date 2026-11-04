#ifndef REAL_TRACK_H
#define REAL_TRACK_H

#include <vector>
#include <cmath>
#include "track.h"

// One physical particle represented by two straight track segments:
// one before the magnet and one after it.
class real_track {
private:
    track trackBefore;
    track trackAfter;
    double momentum;
    double mass;
    double beta;

public:
    // Momentum is obtained during track stitching.
    // Mass is assigned later by particle identification.
    real_track(const track& before, const track& after, double p)
        : trackBefore(before),
          trackAfter(after),
          momentum(p),
          mass(NAN),
          beta(NAN)
    {
        trackBefore.setMom(momentum);
        trackAfter.setMom(momentum);
    }

    void setMom(double p)
    {
        momentum = p;
        trackBefore.setMom(p);
        trackAfter.setMom(p);
    }

    double getMom() const { return momentum; }

    void setMass(double value) { mass = value; }
    double getMass() const { return mass; }

    // Call after setMass().
    void calculateBeta()
    {
        if (!std::isfinite(momentum) || !std::isfinite(mass)){
            beta = NAN;
        } else {
            beta = momentum / std::sqrt(momentum * momentum + mass * mass);
        }

        trackBefore.setBeta(beta);
        trackAfter.setBeta(beta);
    }

    double getBeta() const { return beta; }

    const track& getTrackBefore() const { return trackBefore; }
    const track& getTrackAfter() const { return trackAfter; }
};

// The matching of the two track collections and momentum calculation
// create a real_track(before, after, momentum).
inline std::vector<real_track> stitchTracks(
    const std::vector<track>& tracksBefore,
    const std::vector<track>& tracksAfter)
{
    (void)tracksBefore;
    (void)tracksAfter;
    return {};
}

#endif // REAL_TRACK_H
