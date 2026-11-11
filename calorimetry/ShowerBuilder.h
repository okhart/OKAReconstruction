#ifndef SHOWER_BUILDER_H
#define SHOWER_BUILDER_H

#include <vector>
#include "Hit.h"
#include "Shower.h"

class ShowerBuilder {
public:
    ShowerBuilder();
    ~ShowerBuilder() = default;

    void SetHits(const std::vector<Hit>& hits);  //Load all hits of a detector

    void BuildShowers();   //Make separate Showers

    const std::vector<Shower>& GetShowers() const;    //return Showers

    void Clear();   //Make ready for the next event

private:
    std::vector<Hit> hits;
    std::vector<Shower> showers;
};

#endif