#ifndef SHOWER_H
#define SHOWER_H

#include <vector>
#include "Hit.h"
#include "ParticleIdentifier.h"

class Shower {
public:
    Shower();
    ~Shower() = default;

    void AddHit(const CalHit& hit);
    void RemoveHit(int index);

    const std::vector<CalHit>& GetHits() const;

    int GetMultiplicity() const;

    double GetEnergy() const;

    double GetX() const;
    double GetY() const;
    double GetZ() const;

    double GetTime() const;

    void Identify();    //Identify particle and write into this.type partice type



private:
    std::vector<CalHit> hits;

    double energy;

    double x;
    double y;
    double z;

    double time;

    int multiplicity;

    int type;   //particle
};

#endif