#include <cmath>
#include "limits.hpp"

void Limits::update(const float value) {
    if (_min > value) {
        _min = value;
    }
    if (_max < value) {
        _max = value;
    }
}

std::ostream& operator << (std::ostream& out, const Limits &limits) {
    out << std::floor(limits._min) << ' ' << std::ceil(limits._max);
    return out;
}
