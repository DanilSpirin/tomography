#pragma once

#include "ray.hpp"
#include "distribution.hpp"


class Integration {
public:
    virtual float operator()(const BaseRay &ray, const ElectronDensityDistribution &model) = 0;
    virtual ~Integration() = default;
};

class Trapezium : public Integration {
public:
    float operator() (const BaseRay &ray, const ElectronDensityDistribution &model) override;
    ~Trapezium() {}
};

class Rectangle : public Integration {
public:
    float operator() (const BaseRay &ray, const ElectronDensityDistribution &model) override;
    ~Rectangle() {}
};
