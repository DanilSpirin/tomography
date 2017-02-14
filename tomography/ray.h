#pragma once

#include "point.h"

class BaseRay {
public:
    double integral, time;
    point station, satellite;
};

class Ray : public BaseRay{
public:
    double thetta, phi, angle;
    Ray(const point station, const point satellite, const float time);
    Ray();
    point cross;
    void computeParameters();
private:
    void computeCross();
    void computeAngles();
};

std::istream& operator >> (std::istream&, Ray&);
std::ostream& operator << (std::ostream&, Ray);