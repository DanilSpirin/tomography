#pragma once

#include "vector_sparse.h"
#include "grid.h"
#include "ray.h"
#include "reconstruction.h"
#include "distribution.h"
#include <vector>
#include <list>

// Формирование матрицы задачи из данных
void dataToSle(const std::vector<std::vector<Ray>> &data, std::vector<VectorSparse> &phi, std::vector<double> &integrals, const Grid &test);
void dataToSle(const std::vector<std::vector<Ray>> &data, std::vector<VectorSparse> &phi, const Grid &test);

// Рассчет невязки системы
double computeResidual(const Grid &x, const std::vector<VectorSparse> &A, const std::vector<double> &m);
std::vector<double> computeVectorResidual(const Grid &x, const std::vector<VectorSparse> &A, const std::vector<double> &m);

// Критерий отбора для файлов
int sel(const struct dirent *d);

// Сканирование указанной директории
std::vector<std::vector<Ray>> getData(const char *p, int timeStart = 0, int timeFinish = 24);

// Координаты станций
std::list<std::pair<double, double>> getStationList(std::vector<std::vector<Ray>> data);


void solveSle(Grid &grid, const std::vector<VectorSparse> &matrix, const std::vector<double> integrals, double error, bool onlyPositive = true);
void computeParametrs(Grid &crude, Grid &accurate, std::vector<VectorSparse> sleMatrix, std::vector<double> integrals, bool useSecondGrid, ElectronDensityDistribution &model, Dimension latitude, Dimension longitude, Dimension time, int intervals, int intervalsTime, double initialResidual);

double degreeToRadian(double degree);
double radianToDegree(double radian);

std::list<unsigned> createListOfIntervals(unsigned first, unsigned last);