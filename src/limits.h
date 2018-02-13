#pragma once

#include <limits>
#include <iostream>

class Limits {
    float _min;
    float _max;
public:
    Limits() : _min(std::numeric_limits<float>::max()), _max(std::numeric_limits<float>::lowest()) {}
    Limits(const float min, const float max) : _min(min), _max(max) {}
    void update(const float value);
    friend std::ostream& operator << (std::ostream& out, const Limits &limits);
};
