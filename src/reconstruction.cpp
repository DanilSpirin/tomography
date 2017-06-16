#include <iostream>
#include <vector>

#include "reconstruction.h"
#include "grid.h"

void iterationArt(Grid &x, const std::vector<VectorSparse> &a, const std::vector<double> &m, bool onlyPositive) {
    for (unsigned i = 0; i < m.size(); ++i) {
        double aa = 0;
        double ax = 0;
        for (int k = 0; k < a[i].size(); ++k) {
            const auto element = a[i][k];
            aa += element.value * element.value;
            ax += element.value * x[element.index];
        }
        if (aa == 0) {
            continue;
        }
        const double t = (m[i] - ax) / aa;
        for (unsigned j = 0; j < a[i].size(); ++j) {
            const auto element = a[i][j];
            x[element.index] += element.value * t;
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
            const auto element = a[i][j];
            ax += element.value * x[element.index];
        }
        for (unsigned j = 0; j < a[i].size(); ++j) {
            const auto element = a[i][j];
            dx[element.index] += (m[i] - ax) * element.value;
        }
    }
    // Вычисление коэффициента t для минимизации невязки
    std::vector<double> adx;
    std::vector<double> axy;
    for (unsigned j = 0; j < m.size(); ++j) {
        double adxj = 0;
        double axyj = 0;
        for (unsigned k = 0; k < a[j].size(); ++k) {
            const auto element = a[j][k];
            adxj += element.value * dx[element.index];
            axyj += element.value * x[element.index];
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
