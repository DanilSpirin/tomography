#include <iostream>
#include <cmath>

#include "integration.h"


double Trapezium::operator()(const BaseRay ray, const ElectronDensityDistribution *model) {
    double r = sqrt((ray.satellite - ray.station)*(ray.satellite - ray.station));
    double m = 100;
    double h = r/m;
    point dR = (ray.satellite - ray.station) * (1/m);
    point R = ray.station;
    double value = 0.0;
    
    for (int i = 0; i < m; ++i){
        if (i == 0 || i == m-1)
            value += (*model)(R, ray.time)*h/2;
        else
            value += (*model)(R, ray.time)*h;
        R = R + dR;
    }
    return value;
}

double Rectangle::operator()(const BaseRay ray, const ElectronDensityDistribution *model) {
    double r = sqrt((ray.satellite - ray.station)*(ray.satellite - ray.station));
    double m = 100;
    double h = r/m;
    point dR = (ray.satellite - ray.station) * (1/m);
    point R = ray.station + dR * 0.5;
    double value = 0.0;
    
    for (int i = 0; i < m; ++i){
        value += (*model)(R, ray.time) * h;
        R = R + dR;
    }
    return value;
}