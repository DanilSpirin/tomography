#include <algorithm>
#include <numeric>

#include "reconstruction.h"


void Art::operator()(Grid &x, const std::vector<VectorSparse> &a, const std::vector<double> &m, bool onlyPositive) const {
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
            std::transform(x.begin(), x.end(), x.begin(),
                [] (double value) { return std::max(value, 0.0); }
            );
        }
    }

}

void Sirt::operator()(Grid &x, const std::vector<VectorSparse> &a, const std::vector<double> &m, bool onlyPositive) const {
    // Calculating x increment
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

    // Calculating t coefficient for residual minimization
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

    // Numerator is a sum of axy and adx product values
    // While denominator is a sum of adx values squared
    const double numerator = std::inner_product(axy.begin(), axy.end(), adx.begin(), 0.0);
    const double denominator = std::inner_product(adx.begin(), adx.end(), adx.begin(), 0.0);
    const double t = - numerator / denominator;

    std::transform(dx.begin(), dx.end(), dx.begin(), [t](auto i) { return t * i;});
    std::transform(x.begin(), x.end(), dx.begin(), x.begin(), std::plus<double>());

    if (onlyPositive) {
        std::transform(x.begin(), x.end(), x.begin(),
            [] (double value) { return std::max(value, 0.0); }
        );
    }
}
