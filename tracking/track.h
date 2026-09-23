#ifndef TRACK_H
#define TRACK_H

#include <TFile.h>
#include <TTree.h>
#include <TVector.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <TH2F.h>

using namespace std;

// =====================================================
// ===== PLANE TYPE ENUM ================================
// =====================================================

/// Specifies which coordinates a hit provides
enum class PlaneType { ZX, ZY, ZXY };

// =====================================================
// ===== HIT CLASS =====================================
// =====================================================

/**
 * Represents a single measurement on a detector plane.
 * Stores position (z, x, y), measurement uncertainty (sigma),
 * and the type of plane that produced the hit.
 */
class hit {
private:
    double z;       ///< Z-coordinate (common to all planes)
    double x;       ///< X-coordinate (for ZX and ZXY planes)
    double y;       ///< Y-coordinate (for ZY and ZXY planes)
    double sigma;   ///< Measurement uncertainty
    PlaneType type; ///< Type of detector plane

public:
    hit() : z(0), x(0), y(0), sigma(1e-15), type(PlaneType::ZX) {}

    // ===== SETTERS =====================================
    void setZ(double z_val) { z = z_val; }
    void setX(double x_val) { x = x_val; }
    void setY(double y_val) { y = y_val; }
    void setSigma(double sigma_val) { sigma = sigma_val; }
    void setType(PlaneType plane_type) { type = plane_type; }

    // ===== GETTERS =====================================
    double getZ() const { return z; }
    double getX() const { return x; }
    double getY() const { return y; }
    double getSigma() const { return sigma; }
    PlaneType getType() const { return type; }

    // ===== COORDINATE CHECKS ===========================
    bool hasX() const {
        return (type == PlaneType::ZX || type == PlaneType::ZXY) && !std::isnan(x);
    }
    bool hasY() const {
        return (type == PlaneType::ZY || type == PlaneType::ZXY) && !std::isnan(y);
    }

    /// Returns true if this hit measures the specified axis
    bool hasCoord(PlaneType axis) const {
        if (axis == PlaneType::ZX) return hasX();
        if (axis == PlaneType::ZY) return hasY();
        return hasX() || hasY();
    }

    /// Returns the coordinate value for the specified axis
    double getCoord(PlaneType axis) const {
        if (axis == PlaneType::ZX) return getX();
        if (axis == PlaneType::ZY) return getY();
        return NAN;
    }
};

// =====================================================
// ===== TRACK CLASS ===================================
// =====================================================

/**
 * Represents a reconstructed particle trajectory.
 * Fits straight lines in both ZX and ZY projections using
 * weighted least squares based on hit uncertainties.
 */
class track {
private:
    double a, b;        ///< ZX line: x = a*z + b
    double c, d;        ///< ZY line: y = c*z + d
    const hit** hits;   ///< Array of hits belonging to this track
    int nHits;          ///< Number of hits
    double chi2;        ///< Chi-squared/ndf goodness-of-fit

    // =====================================================
    // ===== WEIGHTED LINE FITTING =========================
    // =====================================================

    /**
     * Performs a weighted linear regression: y = k*x + m.
     * Weights are 1/sigma^2 (inverse variance).
     */
    void fitLineWeighted(const double* x_vals, const double* y_vals,
                         const double* sigma_vals, int n, double& k, double& m) {
        if (n < 2) {
            k = 0;
            m = 0;
            return;
        }

        double sum_w = 0, sum_wx = 0, sum_wy = 0, sum_wxx = 0, sum_wxy = 0;

        for (int i = 0; i < n; i++) {
            double w = 1.0 / (sigma_vals[i] * sigma_vals[i]);
            sum_w  += w;
            sum_wx += w * x_vals[i];
            sum_wy += w * y_vals[i];
            sum_wxx += w * x_vals[i] * x_vals[i];
            sum_wxy += w * x_vals[i] * y_vals[i];
        }

        double denominator = sum_w * sum_wxx - sum_wx * sum_wx;
        if (fabs(denominator) < 1e-10) {
            k = 0;
            m = sum_wy / sum_w;
            return;
        }

        k = (sum_w * sum_wxy - sum_wx * sum_wy) / denominator;
        m = (sum_wy * sum_wxx - sum_wx * sum_wxy) / denominator;
    }

    // =====================================================
    // ===== CHI-SQUARED CALCULATION =======================
    // =====================================================

    /// Computes the total χ² by summing residuals for all measured coordinates.
    void calculateChi2() {
        chi2 = 0;
        for (int i = 0; i < nHits; i++) {
            double sigma = hits[i]->getSigma();

            if (hits[i]->hasX()) {
                double residual = hits[i]->getX() - predictX(hits[i]->getZ());
                chi2 += (residual * residual) / (sigma * sigma);
            }
            if (hits[i]->hasY()) {
                double residual = hits[i]->getY() - predictY(hits[i]->getZ());
                chi2 += (residual * residual) / (sigma * sigma);
            }
        }
        chi2 = (double)chi2/(nHits+1-4);    //Match of 2 tracks in different
                                            //planes (summary, nHits+1 hits
                                            //and 4 parameters)
    }

public:
    // ===== PREDICTION METHODS ============================
    double getA() const { return a; }
    double getB() const { return b; }
    double getC() const { return c; }
    double getD() const { return d; }
    double getChi2() const { return chi2; }

    double predictX(double z) const { return a * z + b; }
    double predictY(double z) const { return c * z + d; }

    // ===== CONSTRUCTORS ==================================
    track() : a(0), b(0), c(0), d(0), hits(nullptr), nHits(0), chi2(0) {}

    track(const hit** hitArray, int n) : hits(nullptr), nHits(n), chi2(0) {
        hits = new const hit*[nHits];
        for (int i = 0; i < nHits; i++) {
            hits[i] = hitArray[i];
        }
        fit();  // Perform the fit immediately
    }

    // ===== DESTRUCTOR ====================================
    ~track() {
        delete[] hits;
    }

    // ===== ACCESSORS =====================================
    int getNHits() const { return nHits; }
    const hit* getHit(int i) const {
        if (i >= 0 && i < nHits) return hits[i];
        return nullptr;
    }

    // =====================================================
    // ===== MAIN FITTING ROUTINE ==========================
    // =====================================================

    /// Fits separate lines for ZX and ZY projections.
    void fit() {
        if (nHits < 2) {
            cerr << "Error: Need at least 2 hits to fit a track" << endl;
            a = b = c = d = 0;
            chi2 = 0;
            return;
        }

        // ---- Fit ZX projection ----
        int zx_count = 0;
        for (int i = 0; i < nHits; i++) {
            if (hits[i]->hasX()) zx_count++;
        }

        if (zx_count >= 2) {
            double* z_vals = new double[zx_count];
            double* x_vals = new double[zx_count];
            double* sigma_vals = new double[zx_count];

            int idx = 0;
            for (int i = 0; i < nHits; i++) {
                if (hits[i]->hasX()) {
                    z_vals[idx] = hits[i]->getZ();
                    x_vals[idx] = hits[i]->getX();
                    sigma_vals[idx] = hits[i]->getSigma();
                    idx++;
                }
            }

            fitLineWeighted(z_vals, x_vals, sigma_vals, zx_count, a, b);

            delete[] z_vals;
            delete[] x_vals;
            delete[] sigma_vals;
        } else {
            a = 0;
            b = 0;
        }

        // ---- Fit ZY projection ----
        int zy_count = 0;
        for (int i = 0; i < nHits; i++) {
            if (hits[i]->hasY()) zy_count++;
        }

        if (zy_count >= 2) {
            double* z_vals = new double[zy_count];
            double* y_vals = new double[zy_count];
            double* sigma_vals = new double[zy_count];

            int idx = 0;
            for (int i = 0; i < nHits; i++) {
                if (hits[i]->hasY()) {
                    z_vals[idx] = hits[i]->getZ();
                    y_vals[idx] = hits[i]->getY();
                    sigma_vals[idx] = hits[i]->getSigma();
                    idx++;
                }
            }

            fitLineWeighted(z_vals, y_vals, sigma_vals, zy_count, c, d);

            delete[] z_vals;
            delete[] y_vals;
            delete[] sigma_vals;
        } else {
            c = 0;
            d = 0;
        }

        calculateChi2();
    }

    // =====================================================
    // ===== PRINT TRACK INFORMATION =======================
    // =====================================================

    void print() const {
        cout << "Track with " << nHits << " hits:" << endl;
        cout << "  ZX: x = " << a << "*z + " << b << endl;
        cout << "  ZY: y = " << c << "*z + " << d << endl;

        int ndf = 0;
        for (int i = 0; i < nHits; i++) {
            if (hits[i]->hasX()) ndf++;
            if (hits[i]->hasY()) ndf++;
        }
        ndf -= 4;  // subtract 4 fit parameters (a,b,c,d)
        if (ndf > 0) {
            cout << "  χ² = " << chi2 * ndf << ", χ²/ndf = " << chi2 << endl;
        } else {
            cout << "  χ² = " << chi2 << endl;
        }

        for (int i = 0; i < nHits; i++) {
            cout << "    Hit " << i << ": Z=" << hits[i]->getZ();
            if (hits[i]->hasX()) cout << ", X=" << hits[i]->getX();
            if (hits[i]->hasY()) cout << ", Y=" << hits[i]->getY();
            cout << ", Type=";
            if (hits[i]->getType() == PlaneType::ZX) cout << "ZX";
            else if (hits[i]->getType() == PlaneType::ZY) cout << "ZY";
            else cout << "ZXY";
            cout << ", Sigma=" << hits[i]->getSigma() << endl;
        }
    }
};

#endif // TRACK_H