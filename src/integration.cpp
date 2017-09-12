#include <cmath>

#include "integration.h"

double Trapezium::operator()(const BaseRay &ray, const ElectronDensityDistribution &model) {
    const double r = sqrt((ray.satellite - ray.station) * (ray.satellite - ray.station));
    const double m = 100;
    const double h = r / m;
    const point dR = (ray.satellite - ray.station) / m;
    point R = ray.station;
    double value = 0.0;

    for (int i = 0; i < m; ++i) {
        if (i == 0 || i == m - 1) {
            value += model(R, ray.time) * h / 2;
        } else {
            value += model(R, ray.time) * h;
        }
        R = R + dR;
    }
    return value;
}

double Rectangle::operator()(const BaseRay &ray, const ElectronDensityDistribution &model) {
    const double r = sqrt((ray.satellite - ray.station) * (ray.satellite - ray.station));
    const double m = 100;
    const double h = r / m;
    const point dR = (ray.satellite - ray.station) / m;
    point R = ray.station + dR * 0.5;
    double value = 0.0;

    for (int i = 0; i < m; ++i) {
        value += model(R, ray.time) * h;
        R = R + dR;
    }
    return value;
}
