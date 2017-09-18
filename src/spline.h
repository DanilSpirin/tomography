#pragma once
#include <vector>

class Spline{
public:
    Spline(const std::vector<double> x, const std::vector<double> y);
    ~Spline() = default;
    double operator() (const double x) const; // возвращает значение в точке x
private:
    //структура, описывающая сплайн на каждом сегменте сетки
    struct SplineTuple {
        double a, b, c, d, x;
    };
    std::vector<SplineTuple> splines;
    void build(const std::vector<double> x, const std::vector<double> y); // построение сплайна по n точкам (x, y)
};
