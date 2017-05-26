#pragma once

#include <iostream>

class Dimension {
public:
    double left, right;
    Dimension();
    Dimension(const double left, const double right, const unsigned itervals, const bool inDegrees = true);
    Dimension(Dimension const &foo);
    void toRadian();
    void toDegrees();
    void expand(const unsigned extra = 16);
    unsigned size() const;
    double length() const;
    std::vector<unsigned> sequence(const double i) const; // TODO: come up withe a better name
    double get(const double x, const double i) const; // TODO: come up with a better name
    friend std::ostream& operator << (std::ostream& out, const Dimension& dimension);
    Dimension& operator = (const Dimension &foo);
private:
    bool inDegrees;
    double step;
    unsigned intervals;
};