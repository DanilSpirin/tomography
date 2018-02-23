#include <cmath>

#include "integration.hpp"

float Trapezium::operator()(const BaseRay &ray, const ElectronDensityDistribution &model) {
    const float r = std::sqrt((ray.satellite - ray.station) * (ray.satellite - ray.station));
    const unsigned m = 100;
    const float h = r / m;
    const point dR = (ray.satellite - ray.station) / m;
    point R = ray.station;
    float value = 0.0;

    for (unsigned i = 0; i < m; ++i) {
        if (i == 0 || i == m - 1) {
            value += model(R, ray.time) * h / 2;
        } else {
            value += model(R, ray.time) * h;
        }
        R = R + dR;
    }
    return value;
}

float Rectangle::operator()(const BaseRay &ray, const ElectronDensityDistribution &model) {
    const float r = std::sqrt((ray.satellite - ray.station) * (ray.satellite - ray.station));
    const unsigned m = 100;
    const float h = r / m;
    const point dR = (ray.satellite - ray.station) / m;
    point R = ray.station + dR * 0.5;
    float value = 0.0;

    for (unsigned i = 0; i < m; ++i) {
        value += model(R, ray.time) * h;
        R = R + dR;
    }
    return value;
}
