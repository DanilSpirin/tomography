#pragma once

#include <vector>

#include "transformation.h"

class Spot {
public:
    Spot(const point &location, const double peak, const double period, const double intensity, const double size);
    double operator () (const point& R, const double time) const;
private:
    point location;
    double peak, period, intensity, size;
};

class Wave {
public:
    Wave(const point &location, const double start, const double period, const double speed);
    double operator () (const point& R, const double time) const;
private:
    point location;
    double start, period, speed;
    double p(const double r, const double v, const double T) const;
    double f(const double t, const double T, const double n) const;
};

class ElectronDensityDistribution {
public:
    virtual ~ElectronDensityDistribution() = default;
    double operator() (point R, const double t) const;
    CoordinateTransformation *coordinateTransformation;
private:
    virtual double value(const point &R, const double t) const = 0;
};

class ChepmanLayer : public ElectronDensityDistribution {
public:
    ChepmanLayer();
    ~ChepmanLayer();
    void add_spot(const point &location, const double peak, const double period, const double intensity, const double size);
    void add_wave(const point &location, const double peak, const double period, const double speed);
private:
    double value(const point &R, const double t) const;
    double nmin, nm, hm, H, d, dt;
    
    std::vector<Spot> spots;
    std::vector<Wave> waves;
};
