#pragma once

#include <limits>
#include <iostream>

class Limits {
    double _min;
    double _max;
public:
    Limits() : _min(std::numeric_limits<double>::max()), _max(std::numeric_limits<double>::lowest()) {}
    Limits(const double min, const double max) : _min(min), _max(max) {}
    void update(const double value);
    friend std::ostream& operator << (std::ostream& out, const Limits &limits);
};