#include <iostream>
#include <fstream>
#include <math.h>
#include "ray.h"


std::istream& operator >> (std::istream& in, Ray &c) {
    in >> c.station >> c.satellite >> c.time >> c.integral;
    return in;
}
std::ostream& operator << (std::ostream &out, Ray c) {
    out << c.station << '\n' << c.satellite << " " << c.time << " " << c.integral;
    return out;
}

Ray::Ray(const point &station, const point &satellite, const float time) {
    this->station = station;
    this->satellite = satellite;
    this->time = time;
}

void Ray::compute_cross() {
    const point dr = satellite - station;
    const float rdr = station * dr;
    const float t = (-rdr + sqrt(rdr * rdr + (dr * dr) * ((Re + h) * (Re + h) - (station * station)))) / (dr * dr);
    cross = station + t * (satellite - station);
}

void Ray::compute_angles() {
    const point diff = (satellite - station);
    const float a = cross * diff;
    const float c = std::sqrt((cross * cross) * (diff * diff));
    const float ac = std::fmin(a / c, 1.0f);
    angle = std::acos(ac);
    thetta = asin(cross.R[2] / (Re + h));
    phi = atan2(cross.R[1], cross.R[0]);
}

void Ray::compute_parameters() {
    compute_cross();
    compute_angles();
}
