#pragma once

#include "point.hpp"

class BaseRay {
public:
    point station, satellite;
    float time, integral;
    BaseRay(const point &station, const point &satellite, float time, float integral);
};

class Ray : public BaseRay{
public:
    float thetta, phi, angle;
    Ray() = delete;
    Ray(const point &station, const point &satellite, float time, float integral);
private:
    void compute_parameters();
};
