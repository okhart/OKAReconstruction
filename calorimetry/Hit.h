#ifndef HIT_H
#define HIT_H

// =====================================================
// ===== PLANE TYPE ====================================
// =====================================================

enum class PlaneType { ZX, ZY, ZXY };

// =====================================================
// ===== DETECTOR ======================================
// =====================================================

enum class Detector { BGD, GAMS, GDA, DT, HODO};

// =====================================================
// ===== HIT CLASS =====================================
// =====================================================
/*
 * Represents one detector measurement.
 * A hit contains z and one or both transverse coordinates.
 */
 class Hit {
    private:
        double z;
        double x;
        double y;
        double sigma;
        PlaneType type;
        Detector detectorID;
        int cellID;
        double time;
    
    public:
        hit() 
            : z(0), x(0), y(0), sigma(1e-15),
            time(0)
        {}
    
        explicit hit(PlaneType type)
            : z(0), x(0), y(0), sigma(1e-15),
              type(type), time(0)
        {}
    
        void setZ(double value) { z = value; }
        void setX(double value) { x = value; }
        void setY(double value) { y = value; }
        void setSigma(double value) { sigma = value; }
        void setType(PlaneType value) { type = value; }
        void setDetector(Detector value) { detectorID=value; }
        void setCellID(int value) { cellID=value; }
        void setTime(double value) { y = value; }
        
    
        double getZ() const { return z; }
        double getX() const { return x; }
        double getY() const { return y; }
        double getSigma() const { return sigma; }
        PlaneType getType() const { return type; }
        int GetCellID() const;
        double GetTime() const;
    
        bool hasX() const {
            return (type == PlaneType::ZX || type == PlaneType::ZXY) && !std::isnan(x);
        }
    
        bool hasY() const {
            return (type == PlaneType::ZY || type == PlaneType::ZXY) && !std::isnan(y);
        }
    
        bool hasCoord(PlaneType axis) const {
            if (axis == PlaneType::ZX) return hasX();
            if (axis == PlaneType::ZY) return hasY();
            return hasX() || hasY();
        }
    
        double getCoord(PlaneType axis) const {
            if (axis == PlaneType::ZX) return getX();
            if (axis == PlaneType::ZY) return getY();
            return NAN;
        }

    protected:
        static void validateDetectorType(Detector detectorID, PlaneType type) {
            bool valid = false;
    
            if (detectorID == DT) && ((type == ZX)||(type == ZY)) {
                valid = true;
            }
            else if (detectorID != DT) && (type == ZXY){
                valid = true
            }
    
            if (!valid) {
                throw std::invalid_argument(
                    "Invalid PlaneType for this detectorID"
                );
            }
        }
    };  

class CalHit : public Hit {
    private:
        double energy;
    
    public:
        CalHit()
            : hit(PlaneType::ZXY),
        {}
    
        double getEnergy() const;
        void setEnergy(double value) { energy=value; } 
    };

    protected:
        static void validateDetectorType(Detector detectorID, PlaneType type) override {
            bool valid = false;
    
            if (type == ZXY) && (detectorID != HODO) && (detectorID != DT){
                valid = true
            }
    
            if (!valid) {
                throw std::invalid_argument(
                    "Invalid PlaneType for hit from calorimeter (CalHit)"
                );
            }
        }

#endif