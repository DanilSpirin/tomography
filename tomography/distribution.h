#pragma once

#include <vector>

#include "transformation.h"

using namespace std;

class Spot {
public:
    Spot();
    Spot(point location, double peak, double period, double intensity, double size);
    Spot(const Spot &a);
    ~Spot();
    double operator () (const point R, const double time) const;
private:
    point location;
    double peak, period, intensity, size;
};

class Wave {
public:
    Wave();
    Wave(point location, double start, double period, double speed);
    Wave(const Wave &a);
    ~Wave();
    double operator () (const point R, const double time) const;
private:
    point location;
    double period, start, speed;
    double p(const double r, const double v, const double T) const;
    double f(const double t, const double T, const double n) const;
};

class ElectronDensityDistribution {
public:
    double operator() (point R, const double t) const;
    CoordinateTransformation *coordinateTransformation;
private:
    virtual double value(const point R, const double t) const = 0;
};

class ChepmanLayer : public ElectronDensityDistribution {
public:
    ChepmanLayer();
    ~ChepmanLayer();
    void addSpot(point location, double peak, double period, double intensity, double size);
    void addWave(point location, double peak, double period, double speed);
private:
    double value(const point R, const double t) const;
    double nmin, nm, hm, H, d, dt;
    
    vector<Spot> spots;
    vector<Wave> waves;
};