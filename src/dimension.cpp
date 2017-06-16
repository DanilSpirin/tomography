#include <cmath>
#include <vector>
#include <numeric>
#include "dimension.h"
#include "ray.h"


Dimension::Dimension(const double left, const double right, const unsigned intervals, const bool isInDegrees) {
    this->left = left;
    this->right = right;
    this->intervals = intervals;
    this->inDegrees = isInDegrees;
    this->step = (right - left) / intervals;
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

void Dimension::expand(const unsigned extra) {
    unsigned m = (unsigned)ceil(extra / step);
    double D = m * step;
    this->left -= D;
    this->right += D;
    this->intervals += 2 * m;
}

unsigned Dimension::size() const {
    return this->intervals;
}

double Dimension::length() const {
    return this->right - this->left;
}

std::vector<unsigned> Dimension::sequence(const double x) const {
    const unsigned i = (unsigned)floor((x - this->left) / this->step);
    const unsigned start_point = i < 2 ? 0 : (i - 1);
    const unsigned end_point = (i + 2) > this->intervals ? this->intervals : (i + 2);
    std::vector<unsigned> result(end_point - start_point + 1);
    std::iota(result.begin(), result.end(), start_point);
    return result;
}

double Dimension::get(const double x, const double i) const {
    return (x - (this->left + i * this->step)) / this->step;
}

std::ostream& operator << (std::ostream& out, const Dimension& dimension) {
    out << dimension.left << ' ' << dimension.right;
    return out;
}