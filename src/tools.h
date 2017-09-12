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

// Сканирование указанной директории
std::vector<std::vector<Ray>> get_data(const std::string &path, const unsigned startTime = 0, const unsigned finishTime = 24);

// Координаты станций
std::list<std::pair<double, double>> getStationList(std::vector<std::vector<Ray>> data);


void solveSle(Grid &grid, const std::vector<VectorSparse> &matrix, const std::vector<double> &integrals, const double error, const bool onlyPositive = true);
void computeParametrs(Grid &crude, Grid &accurate, const std::vector<VectorSparse> &sleMatrix, const std::vector<double> &integrals, const bool useSecondGrid, ElectronDensityDistribution &model, Dimension latitude, Dimension longitude, Dimension time, unsigned intervals, unsigned intervalsTime, double initialResidual);

double degreeToRadian(const double degree);
double radianToDegree(const double radian);

std::list<unsigned> createListOfIntervals(const unsigned first, const unsigned last);
