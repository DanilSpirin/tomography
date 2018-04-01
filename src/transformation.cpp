#include "transformation.hpp"

#include <cmath>

void DecartToGeographic::forward(point &R) {
    point r;
    r.R[2] = std::sqrt(R.R[0] * R.R[0] + R.R[1] * R.R[1] + R.R[2] * R.R[2]);
    r.R[0] = std::atan2(R.R[1], R.R[0]);
    r.R[1] = std::asin(R.R[2] / r.R[2]);
    R = r;
}

void DecartToGeographic::backward(point &R) {
    point r;
    r.R[0] = R.R[2] * std::cos(R.R[1]) * std::cos(R.R[0]);
    r.R[1] = R.R[2] * std::cos(R.R[1]) * std::sin(R.R[0]);
    r.R[2] = R.R[2] * std::sin(R.R[1]);
    R = r;
}
