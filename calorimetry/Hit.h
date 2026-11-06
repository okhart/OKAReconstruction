#ifndef HIT_H
#define HIT_H

class Hit {
public:
    Hit();
    ~Hit() = default;

    int GetDetectorID() const;  //1 = BGD, 2=GAMS, 3 = GDA 
    int GetCellID() const;

    double GetEnergy() const;
    double GetTime() const;

    double GetX() const;
    double GetY() const;
    double GetZ() const;

private:
    int detectorID;
    int cellID;

    double energy;
    double time;

    double x;
    double y;
    double z;
};

#endif