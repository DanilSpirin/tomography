#pragma once

#include <vector>

class Spline {
public:
    Spline(const std::vector<float> x, const std::vector<float> y);
    ~Spline() = default;
    float operator() (const float x) const; // return value at point x
private:
    // structure that defines spline on each segment
    struct SplineTuple {
        float a, b, c, d, x;
    };
    std::vector<SplineTuple> splines;
};
