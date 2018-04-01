#include "point.hpp"

#include <cmath>

point::point(const float x, const float y, const float z) : R {x, y, z} {}

float point::length() const {
    return std::sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2]);
}

float point::length_squared() const {
    return R[0] * R[0] + R[1] * R[1] + R[2] * R[2];
}

float point::radius_squared() const {
    return R[0] * R[0] + R[1] * R[1];
}

point operator + (const point& a, const point& b) {
    return point(a.R[0] + b.R[0], a.R[1] + b.R[1], a.R[2] + b.R[2]);
}

point operator - (const point& a, const point& b) {
    return point(a.R[0] - b.R[0], a.R[1] - b.R[1], a.R[2] - b.R[2]);
}

float operator * (const point& a, const point& b) {
    return a.R[0] * b.R[0] + a.R[1] * b.R[1] + a.R[2] * b.R[2];
}

point operator * (const point& a, const float& b) {
    return point(a.R[0] * b, a.R[1] * b, a.R[2] * b);
}

point operator * (const float& a, const point& b) {
    return point(a * b.R[0], a * b.R[1], a * b.R[2]);
}

point operator / (const point& a, const float& b) {
    return point(a.R[0] / b, a.R[1] / b, a.R[2] / b);
}

std::ostream& operator << (std::ostream& out, point& a) {
    out << a.R[0] << " " << a.R[1] << " " << a.R[2];
    return out;
}

std::istream& operator >> (std::istream& in, point& a) {
    in >> a.R[0] >> a.R[1] >> a.R[2];
    return in;
}
