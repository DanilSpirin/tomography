#include <iostream>
#include <fstream>
#include <cmath>
#include "ray.hpp"


BaseRay::BaseRay(const point &station, const point &satellite, float time, float integral)
        : station(station), satellite(satellite), time(time), integral(integral) {}

void Ray::compute_parameters() {
    const point dr = satellite - station;
    const float rdr = station * dr;
    const float t = (-rdr + std::sqrt(rdr * rdr + (dr * dr) * ((Re + h) * (Re + h) - (station * station)))) / (dr * dr);
    const point cross = station + t * dr;

    const float a = cross * dr;
    const float c = std::sqrt((cross * cross) * (dr * dr));
    const float ac = std::fmin(a / c, 1.0f);
    angle = std::acos(ac);
    thetta = std::asin(cross.R[2] / (Re + h));
    phi = std::atan2(cross.R[1], cross.R[0]);
}

Ray::Ray(const point &station, const point &satellite, const float time, const float integral)
        : BaseRay(station, satellite, time, integral) {
    this->compute_parameters();
}
