#pragma once

#include <iostream>
static const double Re = 6371; // Mean Earth radius
static const double pi = 3.14159265; // PI
static const double h = 300; // Single layer height above Earth

class point {
public:
    double R[3];
    point(const double x = 0, const double y = 0, const double z = Re);
    double length() const;
    double length_squared() const;
    double radius_squared() const;
};

point operator - (const point&a, const point& b);
point operator + (const point&a, const point& b);
point operator * (const double& a, const point& b);
point operator * (const point& a, const double& b);
point operator / (const point& a, const double& b);
double operator * (const point &a, const point& b);

std::ostream& operator << (std::ostream& out, point& a);
std::istream& operator >> (std::istream& in, point& a);
