#pragma once

#include "vector_sparse.h"
#include "grid.h"
#include "ray.h"
#include "reconstruction.h"
#include "distribution.h"
#include <vector>
#include <list>

using namespace std;


// Формирование матрицы задачи из данных
void dataToSle(const vector<vector<Ray>> &data, vector<VectorSparse> &phi, vector<double> &integrals, const Grid &test);
void dataToSle(const vector<vector<Ray>> &data, vector<VectorSparse> &phi, const Grid &test);

// Рассчет невязки системы
double computeResidual(const Grid &x, const vector<VectorSparse> &A, const vector<double> &m);
vector<double> computeVectorResidual(const Grid &x, const vector<VectorSparse> &A, const vector<double> &m);

// Критерий отбора для файлов
int sel(const struct dirent *d);

// Сканирование указанной директории
vector<vector<Ray>> getData(const char *p, int timeStart = 0, int timeFinish = 24);

// Координаты станций
list<pair<double, double>> getStationList(vector<vector<Ray>> data);


void solveSle(Grid &grid, const vector<VectorSparse> &matrix, const vector<double> integrals, double error, bool onlyPositive = true);
void computeParametrs(Grid &crude, Grid &accurate, vector<VectorSparse> sleMatrix, vector<double> integrals, bool useSecondGrid, ElectronDensityDistribution &model, Dimension latitude, Dimension longitude, Dimension time, int intervals, int intervalsTime, double initialResidual);

double degreeToRadian(double degree);
double radianToDegree(double radian);

list<int> createListOfIntervals(int first, int last);