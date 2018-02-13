#pragma once

#include "point.h"

class BaseRay {
public:
    float integral, time;
    point station, satellite;
};

class Ray : public BaseRay{
public:
    float thetta, phi, angle;
    Ray(const point &station, const point &satellite, const float time);
    Ray() = default;
    point cross;
    void compute_parameters();
private:
    void compute_cross();
    void compute_angles();
};

std::istream& operator >> (std::istream&, Ray&);
std::ostream& operator << (std::ostream&, Ray);
