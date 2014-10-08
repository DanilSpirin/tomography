//
//  spline.cpp
//  vector_3D_task
//
//  Created by Даниил Спирин on 11.05.12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//
#include <iostream>
#include "spline.h"

Spline::Spline() : splines (NULL) {};

Spline::~Spline() {
    freeMemory();
}

void Spline::freeMemory() {
    delete[] splines;
    splines = NULL;
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
    splines[0].c = splines[n-1].c = 0;
    
    // Метод прогонки
    // Нахождение прогоночных коэффициентов
    double *alpha = new double[n - 1];
    double *beta = new double[n - 1];
    double A, B, C, F, h_i, h_i1, z;
    alpha[0] = beta[0] = 0.;
    for (std::size_t i = 1; i < n - 1; ++i) {
        h_i = x[i] - x[i - 1], h_i1 = x[i + 1] - x[i];
        A = h_i;
        C = 2. * (h_i + h_i1);
        B = h_i1;
        F = 6. * ((y[i + 1] - y[i]) / h_i1 - (y[i] - y[i - 1]) / h_i);
        z = (A * alpha[i - 1] + C);
        alpha[i] = -B / z;
        beta[i] = (F - A * beta[i - 1]) / z;
    }
    
    // Нахождение решения - обратный ход прогонки
    for (std::size_t i = n - 2; i > 0; --i)
        splines[i].c = alpha[i] * splines[i + 1].c + beta[i];
    
    delete[] beta;
    delete[] alpha;
    
    // Нахождение оставшихся коэффициентов
    for (std::size_t i = n - 1; i > 0; --i) {
        double h_i = x[i] - x[i - 1];
        splines[i].d = (splines[i].c - splines[i - 1].c) / h_i;
        splines[i].b = h_i * (2. * splines[i].c + splines[i - 1].c) / 6. + (y[i] - y[i - 1]) / h_i;
    }
}

double Spline::operator()(double x) {
    if (!splines)
        return std::numeric_limits<double>::quiet_NaN();
    
    SplineTuple *s;
    if (x <= splines[0].x) // Если х меньше левой границы то возвращаем значение на левой границе
        return splines[0].a;
    else if (x >= splines[n - 1].x) // Аналогично с правой
        return splines[n-1].a;
    else { // Иначе проводим бинарный поиск нужного элемента массива
        std::size_t i = 0, j = n - 1;
        while (i + 1 < j) {
            std::size_t k = i + (j - i) / 2;
            if (x <= splines[k].x)
                j = k;
            else
                i = k;
        }
        s = splines + j;
    }
    
    double dx = (x - s->x);

    return s->a + (s->b + (s->c / 2. + s->d * dx / 6.) * dx) * dx;
}