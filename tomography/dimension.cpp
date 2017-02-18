#include <cmath>
#include <vector>
#include <numeric>
#include "dimension.h"
#include "ray.h"


Dimension::Dimension(const double left, const double right, const size_t intervals, const bool isInDegrees) {
    this->left = left;
    this->right = right;
    this->intervals = intervals;
    this->inDegrees = isInDegrees;
    this->step = (right - left) / intervals;
}

Dimension::Dimension(Dimension const &foo) {
    this->left = foo.left;
    this->right = foo.right;
    this->intervals = foo.intervals;
    this->step = foo.step;
    this->inDegrees = foo.inDegrees;
}

Dimension::Dimension() {

}

void Dimension::toRadian() {
    if (inDegrees) {
        left = left / 180 * pi;
        right = right / 180 * pi;
        step = step / 180 * pi;
        inDegrees = false;
    }
}

void Dimension::toDegrees() {
    if (!inDegrees) {
        left = left / pi * 180;
        right = right / pi * 180;
        step = step / pi * 180;
        inDegrees = true;
    }
}

void Dimension::expand(const uint16_t extra) {
    uint16_t m = (uint16_t)ceil(extra / step);
    double D = m * step;
    this->left -= D;
    this->right += D;
    this->intervals += 2 * m;
}

size_t Dimension::size() const {
    return this->intervals;
}

std::vector<size_t> Dimension::sequence(const double x) const {
    const auto i = (size_t)floor((x - this->left) / this->step);
    const size_t i1 = i < 2 ? 0 : (i - 1);
    const size_t i2 = (i + 2) > this->intervals ? this->intervals : (i + 2);
    std::vector<size_t> result(i2 - i1 + 1);
    std::iota(result.begin(), result.end(), i1);
    return result;
}

double Dimension::get(const double x, const double i) const {
    return (x - (this->left + i * this->step)) / this->step;
}

std::ostream& operator << (std::ostream& out, const Dimension& dimension) {
    out << dimension.left << ' ' << dimension.right;
    return out;
}