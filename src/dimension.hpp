#pragma once

#include <iostream>
#include <vector>

class Dimension {
public:
    float left, right;
    Dimension() = default;
    Dimension(const float left, const float right, const unsigned itervals, const bool inDegrees = true);
    void to_radian();
    void to_degrees();
    void expand(const unsigned extra = 16);
    unsigned size() const;
    std::vector<unsigned> sequence(const float i) const; // TODO: come up withe a better name
    float get(const float x, const float i) const; // TODO: come up with a better name
    friend std::ostream& operator << (std::ostream& out, const Dimension& dimension);
private:
    bool inDegrees;
    float step;
    unsigned intervals;
};
