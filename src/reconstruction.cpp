#include <iostream>
#include <vector>

#include "reconstruction.h"
#include "grid.h"

void iterationArt(Grid &x, const std::vector<VectorSparse> &a, const std::vector<double> &m, bool onlyPositive) {
    for (unsigned i = 0; i < m.size(); ++i) {
        double aa = 0;
        double ax = 0;
        for (unsigned k = 0; k < a[i].size(); ++k) {
            const auto [index, value] = a[i][k];
            aa += value * value;
            ax += value * x[index];
        }
        if (aa == 0) {
            continue;
        }
        const double t = (m[i] - ax) / aa;
        for (unsigned j = 0; j < a[i].size(); ++j) {
            const auto [index, value] = a[i][j];
            x[index] += value * t;
        }
        if (onlyPositive) {
            for (unsigned j = 0; j < x.size(); ++j) {
                if (x[j] < 0) {
                    x[j] = 0;
                }
            }
        }
    }
}

void iterationSirt(Grid &x, const std::vector<VectorSparse> &a, const std::vector<double> &m, bool onlyPositive) {
    // Вычисление приращения dx
    std::vector<double> dx(x.size(), 0);
    for (unsigned i = 0; i < m.size(); ++i) {
        double ax = 0;
        for (unsigned j = 0; j < a[i].size(); ++j) {
            const auto [index, value] = a[i][j];
            ax += value * x[index];
        }
        for (unsigned j = 0; j < a[i].size(); ++j) {
            const auto [index, value] = a[i][j];
            dx[index] += (m[i] - ax) * value;
        }
    }
    // Вычисление коэффициента t для минимизации невязки
    std::vector<double> adx;
    std::vector<double> axy;
    for (unsigned j = 0; j < m.size(); ++j) {
        double adxj = 0;
        double axyj = 0;
        for (unsigned k = 0; k < a[j].size(); ++k) {
            const auto [index, value] = a[j][k];
            adxj += value * dx[index];
            axyj += value * x[index];
        }
        adx.push_back(adxj);
        axy.push_back(axyj - m[j]);
    }
    double sum1 = 0, sum2 = 0;
    for (unsigned j = 0; j < adx.size(); ++j) {
        sum1 += axy[j] * adx[j];
        sum2 += adx[j] * adx[j];
    }
    const double t = -sum1 / sum2;

    // Вычисляем новый x
    for (unsigned k = 0; k < x.size(); ++k) {
        x[k] += t * dx[k];
    }
    if (onlyPositive) {
        for (unsigned k = 0; k < x.size(); ++k) {
            if (x[k] < 0) {
                x[k] = 0;
            }
        }
    }
}
