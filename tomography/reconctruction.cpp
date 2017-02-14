#include <iostream>
#include <vector>

#include "reconstruction.h"
#include "grid.h"

using namespace std;

void iterationArt(Grid &x, const vector<VectorSparse> &a, const vector<double> &m, bool onlyPositive) {
    for (int i = 0; i < m.size(); ++i) {
        double t = 0, aa = 0, ax = 0;
        for (int k = 0; k < a[i].getSize(); ++k) {
            aa += a[i].getPhi(k) * a[i].getPhi(k);
            ax += a[i].getPhi(k) * x[a[i].getNumber(k)];
        }
        if (aa == 0) {
            continue;
        }
        t = (m[i] - ax) / aa;
        for (int j = 0; j < a[i].getSize(); ++j) {
            x[a[i].getNumber(j)] += (a[i].getPhi(j) * t);
        }
        if (onlyPositive) {
            for (int j = 0; j < x.size(); ++j) {
                if (x[j] < 0) {
                    x[j] = 0;
                }
            }
        }
    }
}

void iterationSirt(Grid &x, const vector<VectorSparse> &a, const vector<double> &m, bool onlyPositive) {
    // Вычисление приращения dx
    vector<float> dx(x.size(), 0);
    double t = 0;
    for (int i = 0; i < m.size(); ++i) {
        double ax = 0;
        for (int j = 0; j < a[i].getSize(); ++j) {
            ax += a[i].getPhi(j) * x[a[i].getNumber(j)];
        }
        for (int j = 0; j < a[i].getSize() ; ++j) {
            dx[a[i].getNumber(j)] += (m[i] - ax)*a[i].getPhi(j);
        }
    }
    // Вычисление коэффициента t для минимизации невязки
    vector<float> adx;
    vector<float> axy;
    for (int j = 0; j < m.size(); ++j) {
        double adxj = 0;
        double axyj = 0;
        for (int k = 0; k < a[j].getSize(); ++k) {
            adxj += a[j].getPhi(k) * dx[a[j].getNumber(k)];
            axyj += a[j].getPhi(k) * x[a[j].getNumber(k)];
        }
        adx.push_back(adxj);
        axy.push_back(axyj - m[j]);
    }
    double sum1=0, sum2=0;
    for (int j = 0; j < adx.size(); ++j) {
        sum1 += axy[j]*adx[j];
        sum2 += adx[j]*adx[j];
    }
    t = -sum1 / sum2;
    
    // Вычисляем новый x
    for (int k = 0; k < x.size(); ++k) {
        x[k] += t * dx[k];
    }
    if (onlyPositive) {
        for (int k = 0; k < x.size(); ++k) {
            if (x[k] < 0) {
                x[k] = 0;
            }
        }
    }
}
