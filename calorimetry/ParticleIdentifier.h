#ifndef PARTICLE_IDENTIFIER_H
#define PARTICLE_IDENTIFIER_H

#include "Cluster.h"

enum class ParticleType: uint8_t {
    Unknown = 0,
    Gamma = 1,
    Positron = 2,
    Neutrino = 4,
    Muon = 5,
    Pi0 = 7,
    Pip = 8, //Pi+
    Pim = 9 //Pi-

};

class ParticleIdentifier {
public:
    ParticleIdentifier();
    ~ParticleIdentifier() = default;

    ParticleType Identify(const Cluster& cluster);

private:
    //Some variables
};

#endif