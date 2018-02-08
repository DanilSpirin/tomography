#include <algorithm>
#include <limits>

#include "spline.h"

Spline::Spline(const std::vector<double> x, const std::vector<double> y) {
    auto n = x.size();
    splines.resize(x.size());
    for (std::size_t i = 0; i < n; ++i) {
        splines[i].x = x[i];
        splines[i].a = y[i];
    }

    // Natural boundary conditions
    splines.front().c = splines.back().c = 0;

    std::vector<double> alpha(n - 1);
    std::vector<double> beta(n - 1);

    // Tridiagonal matrix algorithm
    for (std::size_t i = 1; i < n - 1; ++i) {
        const double A = x[i] - x[i - 1];
        const double B = x[i + 1] - x[i];
        const double C = 2.0 * (A + B);
        const double F = 6.0 * ((y[i + 1] - y[i]) / B - (y[i] - y[i - 1]) / A);
        const double z = A * alpha[i - 1] + C;
        alpha[i] = -B / z;
        beta[i] = (F - A * beta[i - 1]) / z;
    }

    for (std::size_t i = n - 2; i > 0; --i) {
        splines[i].c = alpha[i] * splines[i + 1].c + beta[i];
    }

    for (std::size_t i = n - 1; i > 0; --i) {
        const double h = x[i] - x[i - 1];
        splines[i].d = (splines[i].c - splines[i - 1].c) / h;
        splines[i].b = h * (2.0 * splines[i].c + splines[i - 1].c) / 6.0 + (y[i] - y[i - 1]) / h;
    }
}

double Spline::operator()(const double x) const {
    if (!splines.size()) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    // Do not interpolate past boundary values
    if (x <= splines.front().x) {
        return splines.front().a;
    }
    if (x >= splines.back().x) {
        return splines.back().a;
    }

    auto it = std::lower_bound(splines.begin(), splines.end(), x, [](SplineTuple left, double right) {
        return left.x < right;
    });

    const auto dx = x - it->x;
    return it->a + (it->b + (it->c / 2.0 + it->d * dx / 6.0) * dx) * dx;
}
