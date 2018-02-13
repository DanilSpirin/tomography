#pragma once

#include <iostream>
static const float Re = 6371; // Mean Earth radius
static const float pi = 3.14159265f; // PI
static const float h = 300; // Single layer height above Earth

class point {
public:
    float R[3];
    point(const float x = 0, const float y = 0, const float z = Re);
    float length() const;
    float length_squared() const;
    float radius_squared() const;
};

point operator - (const point&a, const point& b);
point operator + (const point&a, const point& b);
point operator * (const float& a, const point& b);
point operator * (const point& a, const float& b);
point operator / (const point& a, const float& b);
float operator * (const point &a, const point& b);

std::ostream& operator << (std::ostream& out, point& a);
std::istream& operator >> (std::istream& in, point& a);
