#include <iostream>
#include <limits>
#include "spline.h"

Spline::Spline() : splines(nullptr) {
}

Spline::~Spline() {
    freeMemory();
}

void Spline::freeMemory() {
    delete[] splines;
    splines = nullptr;
}

void Spline::build(const double *x, const double *y, std::size_t n) {
    freeMemory();
    this->n = n;

    // Задание начальных условий
    splines = new SplineTuple[n];
    for (std::size_t i = 0; i < n; ++i) {
        splines[i].x = x[i];
        splines[i].a = y[i];
    }

    // Граничные условия
    splines[0].c = splines[n - 1].c = 0;

    // Метод прогонки
    // Нахождение прогоночных коэффициентов
    double *alpha = new double[n - 1];
    double *beta = new double[n - 1];
    alpha[0] = beta[0] = 0.;
    for (std::size_t i = 1; i < n - 1; ++i) {
        const double A = x[i] - x[i - 1];
        const double B = x[i + 1] - x[i];
        const double C = 2.0 * (A + B);
        const double F = 6.0 * ((y[i + 1] - y[i]) / B - (y[i] - y[i - 1]) / A);
        const double z = A * alpha[i - 1] + C;
        alpha[i] = -B / z;
        beta[i] = (F - A * beta[i - 1]) / z;
    }

    // Нахождение решения - обратный ход прогонки
    for (std::size_t i = n - 2; i > 0; --i) {
        splines[i].c = alpha[i] * splines[i + 1].c + beta[i];
    }

    delete[] beta;
    delete[] alpha;

    // Нахождение оставшихся коэффициентов
    for (std::size_t i = n - 1; i > 0; --i) {
        const double h = x[i] - x[i - 1];
        splines[i].d = (splines[i].c - splines[i - 1].c) / h;
        splines[i].b = h * (2.0 * splines[i].c + splines[i - 1].c) / 6.0 + (y[i] - y[i - 1]) / h;
    }
}

double Spline::operator()(const double x) const {
    if (!splines)
        return std::numeric_limits<double>::quiet_NaN();

    SplineTuple *s;
    if (x <= splines[0].x) { // Если х меньше левой границы то возвращаем значение на левой границе
        return splines[0].a;
    } else if (x >= splines[n - 1].x) { // Аналогично с правой
        return splines[n - 1].a;
    } else { // Иначе проводим бинарный поиск нужного элемента массива
        std::size_t i = 0, j = n - 1;
        while (i + 1 < j) {
            std::size_t k = i + (j - i) / 2;
            if (x <= splines[k].x) {
                j = k;
            } else {
                i = k;
            }
        }
        s = splines + j;
    }

    const double dx = (x - s->x);

    return s->a + (s->b + (s->c / 2.0 + s->d * dx / 6.0) * dx) * dx;
}