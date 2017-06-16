#pragma once

#include "ray.h"
#include "distribution.h"


class Integration {
public:
    virtual double operator()(const BaseRay &ray, const ElectronDensityDistribution &model) = 0;
};

class Trapezium : public Integration {
public:
    double operator() (const BaseRay &ray, const ElectronDensityDistribution &model);
};

class Rectangle : public Integration {
public:
    double operator() (const BaseRay &ray, const ElectronDensityDistribution &model);
};