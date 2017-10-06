#pragma once
#include <vector>

class Spline{
public:
    Spline(const std::vector<double> x, const std::vector<double> y);
    ~Spline() = default;
    double operator() (const double x) const; // return value at point x
private:
    // structure that defines spline on each segment
    struct SplineTuple {
        double a, b, c, d, x;
    };
    std::vector<SplineTuple> splines;
    void build(const std::vector<double> x, const std::vector<double> y);
};
