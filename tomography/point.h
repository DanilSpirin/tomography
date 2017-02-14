#pragma once

#include <iostream>
static const double Re = 6371; // средний радиус земли
static const double pi = 3.14159265; // число пи
static const double h = 300; // высота тонкого слоя

class point {
public:
    double R[3];
    point(double x = 0, double y = 0, double z = 0);
    point(point const &a);
};

point operator - (const point&a, const point& b);
point operator + (const point&a, const point& b);
point operator * (const double& a, const point& b);
point operator * (const point& a, const double& b);
point operator / (const point& a, const double& b);
double operator * (const point &a, const point& b);

std::ostream& operator << (std::ostream& out, point& a);
std::istream& operator >> (std::istream& in, point& a);