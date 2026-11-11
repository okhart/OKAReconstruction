#ifndef TRACK_H
#define TRACK_H

#include <TFile.h>
#include <TTree.h>
#include <TVector.h>
#include <TH2F.h>
#include <Hit.h>

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

#include "q.h"

using namespace std;

// =====================================================
// ===== TRACK CLASS ===================================
// =====================================================
/**
 * One reconstructed straight segment.
 * The track stores non-owning pointers to the hits used in the fit.
 * The hit objects must outlive the track object.
 */
class track {
private:
    double a, b;                  // x = a*z + b
    double c, d;                  // y = c*z + d
    vector<const Hit*> hits;      // non-owning references to hits
    double chi2;                  // chi2/ndf when ndf > 0
    double momentum;              // assigned after stitching
    double beta;                  // v/c, assigned after PID
    CovarianceByZ covarianceByZ;

    void fitLineWeighted(const double* xValues,
                         const double* yValues,
                         const double* sigmaValues,
                         int n,
                         double& slope,
                         double& intercept) {
        if (n < 2) {
            slope = 0;
            intercept = 0;
            return;
        }

        double sumW = 0;
        double sumWX = 0;
        double sumWY = 0;
        double sumWXX = 0;
        double sumWXY = 0;

        for (int i = 0; i < n; ++i) {
            const double w = 1.0 / (sigmaValues[i] * sigmaValues[i]);
            sumW += w;
            sumWX += w * xValues[i];
            sumWY += w * yValues[i];
            sumWXX += w * xValues[i] * xValues[i];
            sumWXY += w * xValues[i] * yValues[i];
        }

        const double denominator = sumW * sumWXX - sumWX * sumWX;
        if (fabs(denominator) < 1e-10) {
            slope = 0;
            intercept = sumWY / sumW;
            return;
        }

        slope = (sumW * sumWXY - sumWX * sumWY) / denominator;
        intercept = (sumWY * sumWXX - sumWX * sumWXY) / denominator;
    }

    void calculateChi2() {
        chi2 = 0;
        int nMeasurements = 0;

        for (const Hit* h : hits) {
            const double sigma = h->getSigma();

            if (h->hasX()) {
                const double residual = h->getX() - predictX(h->getZ());
                chi2 += (residual * residual) / (sigma * sigma);
                ++nMeasurements;
            }

            if (h->hasY()) {
                const double residual = h->getY() - predictY(h->getZ());
                chi2 += (residual * residual) / (sigma * sigma);
                ++nMeasurements;
            }
        }

        // Four fitted parameters: a, b, c, d.
        const int ndf = nMeasurements - 4;
        if (ndf > 0)
            chi2 /= static_cast<double>(ndf);
    }

public:
    track()
        : a(0), b(0), c(0), d(0), chi2(0),
          momentum(NAN), beta(NAN) {}

    explicit track(const vector<const Hit*>& trackHits)
        : a(0), b(0), c(0), d(0), hits(trackHits), chi2(0),
          momentum(NAN), beta(NAN) {
        fit();
    }

    // Straight-line parameters and fit quality.
    double getA() const { return a; }
    double getB() const { return b; }
    double getC() const { return c; }
    double getD() const { return d; }
    double getChi2() const { return chi2; }

    double predictX(double z) const { return a * z + b; }
    double predictY(double z) const { return c * z + d; }

    // Kinematic parameters are assigned after stitching / PID.
    void setMom(double value) { momentum = value; }
    double getMom() const { return momentum; }

    void setBeta(double value) { beta = value; }
    double getBeta() const { return beta; }

    int getNHits() const { return static_cast<int>(hits.size()); }

    const Hit* getHit(int i) const {
        if (i >= 0 && static_cast<size_t>(i) < hits.size())
            return hits[static_cast<size_t>(i)];
        return nullptr;
    }

    const vector<const Hit*>& getHits() const { return hits; }

    // =====================================================
    // ===== MAIN FITTING ROUTINE ==========================
    // =====================================================
    void fit() {
        if (hits.size() < 2) {
            cerr << "Error: Need at least 2 hits to fit a track" << endl;
            a = b = c = d = 0;
            chi2 = 0;
            return;
        }

        // ---- Fit ZX projection ----
        vector<double> zX, xValues, sigmaX;
        for (const Hit* h : hits) {
            if (!h->hasX()) continue;
            zX.push_back(h->getZ());
            xValues.push_back(h->getX());
            sigmaX.push_back(h->getSigma());
        }

        if (zX.size() >= 2) {
            fitLineWeighted(zX.data(), xValues.data(), sigmaX.data(),
                            static_cast<int>(zX.size()), a, b);
        } else {
            a = 0;
            b = 0;
        }

        // ---- Fit ZY projection ----
        vector<double> zY, yValues, sigmaY;
        for (const Hit* h : hits) {
            if (!h->hasY()) continue;
            zY.push_back(h->getZ());
            yValues.push_back(h->getY());
            sigmaY.push_back(h->getSigma());
        }

        if (zY.size() >= 2) {
            fitLineWeighted(zY.data(), yValues.data(), sigmaY.data(),
                            static_cast<int>(zY.size()), c, d);
        } else {
            c = 0;
            d = 0;
        }

        calculateChi2();
    }

    void print() const {
        cout << "Track with " << hits.size() << " hits:" << endl;
        cout << "  ZX: x = " << a << "*z + " << b << endl;
        cout << "  ZY: y = " << c << "*z + " << d << endl;

        int nMeasurements = 0;
        for (const Hit* h : hits) {
            if (h->hasX()) ++nMeasurements;
            if (h->hasY()) ++nMeasurements;
        }
        const int ndf = nMeasurements - 4;

        if (ndf > 0)
            cout << "  chi2 = " << chi2 * ndf << ", chi2/ndf = " << chi2 << endl;
        else
            cout << "  chi2 = " << chi2 << " (ndf <= 0)" << endl;

        for (size_t i = 0; i < hits.size(); ++i) {
            const Hit* h = hits[i];
            cout << "    Hit " << i << ": Z=" << h->getZ();
            if (h->hasX()) cout << ", X=" << h->getX();
            if (h->hasY()) cout << ", Y=" << h->getY();
            cout << ", Type=";
            if (h->getType() == PlaneType::ZX) cout << "ZX";
            else if (h->getType() == PlaneType::ZY) cout << "ZY";
            else cout << "ZXY";
            cout << ", Sigma=" << h->getSigma() << endl;
        }
    }

    void calculateCovariances(const QByZ& qByZ) {
        // TODO:
        // Run Kalman along this straight track.
        // At each z obtain Q(z) from qByZ.
        // Store resulting covariance in covarianceByZ[z].
        (void)qByZ;
    }

    double getSigmaX(double z) const {
        // TODO: return sqrt(covarianceByZ[z].varX)
        (void)z;
        return NAN;
    }

    double getSigmaY(double z) const {
        // TODO: return sqrt(covarianceByZ[z].varY)
        (void)z;
        return NAN;
    }
};

#endif // TRACK_H
