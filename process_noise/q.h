#ifndef Q_H
#define Q_H

#include <functional>
#include <map>

// TODO: replace with actual Q matrix
struct QMatrix
{
};

using QFunction = std::function<QMatrix(
    double momentum,
    double velocityBeta,
    double alpha,
    double beta
)>;

// z -> Q-function for this detector/material region
using QByZ = std::map<double, QFunction>;

struct TrackCovariance
{
    double varX = NAN;
    double varY = NAN;

    // TODO: add slope/cross-covariances
};

using CovarianceByZ = std::map<double, TrackCovariance>;

#endif