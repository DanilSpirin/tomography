#include <cmath>
#include "limits.h"

void Limits::update(const double value) {
    if (_min > value) {
        _min = value;
    }
    if (_max < value) {
        _max = value;
    }
}

std::ostream& operator << (std::ostream& out, const Limits &limits) {
    out << floor(limits._min) << ' ' << ceil(limits._max);
    return out;
}
