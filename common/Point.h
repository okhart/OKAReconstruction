#ifndef Point_h
#define Point_h
class Point
{
public:
    Point()
    {
        fx = 0;
        fy = 0;
        fz = 0;
    }
    Point(float xx, float yy, float zz)
    {
        fx = xx;
        fy = yy;
        fz = zz;
    }
    float X() { return fx; }
    float Y() { return fy; }
    float Z() { return fz; }
    float x() { return fx; }
    float y() { return fy; }
    float z() { return fz; }
    float GetX() { return fx; }
    float GetY() { return fy; }
    float GetZ() { return fz; }

private:
    float fx;
    float fy;
    float fz;
};
#endif