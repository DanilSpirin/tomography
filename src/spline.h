#pragma once

class Spline{
public:
    Spline();
    ~Spline();
    void build(const double *x, const double *y, std::size_t n); // построение сплайна по n точкам (x, y)
    double operator() (const double x) const; // возвращает значение в точке x
private:
    //структура, описывающая сплайн на каждом сегменте сетки
    struct SplineTuple{
        double a, b, c, d, x;
    };
    SplineTuple *splines;
    std::size_t n;
    void freeMemory();
};