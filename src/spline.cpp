#include <algorithm>
#include <limits>

#include "spline.h"

Spline::Spline(const std::vector<float> x, const std::vector<float> y) {
    auto n = x.size();
    splines.resize(x.size());
    for (std::size_t i = 0; i < n; ++i) {
        splines[i].x = x[i];
        splines[i].a = y[i];
    }

    // Natural boundary conditions
    splines.front().c = splines.back().c = 0;

    std::vector<float> alpha(n - 1);
    std::vector<float> beta(n - 1);

    // Tridiagonal matrix algorithm
    for (std::size_t i = 1; i < n - 1; ++i) {
        const float A = x[i] - x[i - 1];
        const float B = x[i + 1] - x[i];
        const float C = 2.f * (A + B);
        const float F = 6.f * ((y[i + 1] - y[i]) / B - (y[i] - y[i - 1]) / A);
        const float z = A * alpha[i - 1] + C;
        alpha[i] = -B / z;
        beta[i] = (F - A * beta[i - 1]) / z;
    }

    for (std::size_t i = n - 2; i > 0; --i) {
        splines[i].c = alpha[i] * splines[i + 1].c + beta[i];
    }

    for (std::size_t i = n - 1; i > 0; --i) {
        const float h = x[i] - x[i - 1];
        splines[i].d = (splines[i].c - splines[i - 1].c) / h;
        splines[i].b = h * (2.f * splines[i].c + splines[i - 1].c) / 6.f + (y[i] - y[i - 1]) / h;
    }
}

float Spline::operator()(const float x) const {
    if (splines.empty()) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    // Do not interpolate past boundary values
    if (x <= splines.front().x) {
        return splines.front().a;
    }
    if (x >= splines.back().x) {
        return splines.back().a;
    }

    auto it = std::lower_bound(splines.begin(), splines.end(), x, [](SplineTuple left, float right) {
        return left.x < right;
    });

    const auto dx = x - it->x;
    return it->a + (it->b + (it->c / 2.f + it->d * dx / 6.f) * dx) * dx;
}
