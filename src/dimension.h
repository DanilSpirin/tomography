#pragma once

#include <iostream>
#include <vector>

class Dimension {
public:
    double left, right;
    Dimension() = default;
    Dimension(const double left, const double right, const unsigned itervals, const bool inDegrees = true);
    void to_radian();
    void to_degrees();
    void expand(const unsigned extra = 16);
    unsigned size() const;
    double length() const;
    std::vector<unsigned> sequence(const double i) const; // TODO: come up withe a better name
    double get(const double x, const double i) const; // TODO: come up with a better name
    friend std::ostream& operator << (std::ostream& out, const Dimension& dimension);
private:
    bool inDegrees;
    double step;
    unsigned intervals;
};
