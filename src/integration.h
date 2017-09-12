#pragma once

#include "ray.h"
#include "distribution.h"


class Integration {
public:
    virtual double operator()(const BaseRay &ray, const ElectronDensityDistribution &model) = 0;
    virtual ~Integration() = default;
};

class Trapezium : public Integration {
public:
    virtual double operator() (const BaseRay &ray, const ElectronDensityDistribution &model);
    ~Trapezium() {}
};

class Rectangle : public Integration {
public:
    virtual double operator() (const BaseRay &ray, const ElectronDensityDistribution &model);
    ~Rectangle() {}
};
