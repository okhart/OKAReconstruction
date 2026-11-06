#ifndef CLUSTER_H
#define CLUSTER_H

#include <vector>
#include "Hit.h"
#include "ParticleIdentifier.h"

class Cluster {
public:
    Cluster();
    ~Cluster() = default;

    void AddHit(const Hit& hit);
    void RemoveHit(int index);

    const std::vector<Hit>& GetHits() const;

    int GetMultiplicity() const;

    double GetEnergy() const;

    double GetX() const;
    double GetY() const;
    double GetZ() const;

    double GetTime() const;

    void Identify();    //Identify particle and write into this.type partice type



private:
    std::vector<Hit> hits;

    double energy;

    double x;
    double y;
    double z;

    double time;

    int multiplicity;

    int type;   //particle
};

#endif