#pragma once

#include <vector>
#include <memory>

#include "transformation.h"

class Spot {
public:
    Spot(const point &location, const float peak, const float period, const float intensity, const float size);
    float operator () (const point& R, const float time) const;
private:
    point location;
    float peak, period, intensity, size;
};

class Wave {
public:
    Wave(const point &location, const float start, const float period, const float speed);
    float operator () (const point& R, const float time) const;
private:
    point location;
    float start, period, speed;
    float p(const float r, const float v, const float T) const;
    float f(const float t, const float T, const float n) const;
};

class ElectronDensityDistribution {
public:
    virtual ~ElectronDensityDistribution() = default;
    float operator() (point R, const float t) const;
    std::unique_ptr<CoordinateTransformation> coordinateTransformation;
private:
    virtual float value(const point &R, const float t) const = 0;
};

class ChepmanLayer : public ElectronDensityDistribution {
public:
    ChepmanLayer();
    ~ChepmanLayer() override;
    void add_spot(const point &location, const float peak, const float period, const float intensity, const float size);
    void add_wave(const point &location, const float peak, const float period, const float speed);
private:
    float value(const point &R, const float t) const override;
    float nmin, nm, hm, H, d, dt;
    
    std::vector<Spot> spots;
    std::vector<Wave> waves;
};
