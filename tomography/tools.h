//
//  tools.h
//  vector_3D_task
//
//  Created by Даниил Спирин on 19.05.12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef vector_3D_task_tools_h
#define vector_3D_task_tools_h

#include "vector_sparse.h"
#include "grid.h"
#include "ray.h"
#include "reconstruction.h"
#include "distribution.h"
#include <vector>
#include <list>

using namespace std;


// Формирование матрицы задачи из данных
void dataToSle(const vector<vector<Ray> > &data, vector<VectorSparse> &phi, vector<float> &integrals, const Grid &test);
void dataToSle(const vector<vector<Ray> > &data, vector<VectorSparse> &phi, const Grid &test);

// Рассчет невязки системы
double computeResidual(const Grid &x, const vector<VectorSparse> &A, const vector<float> &m);
vector<float> computeVectorResidual(const Grid &x, const vector<VectorSparse> &A, const vector<float> &m);

// Критерий отбора для файлов
int sel(const struct dirent *d);

// Сканирование указанной директории
vector<vector<Ray> > getData(const char *p, int timeStart = 0, int timeFinish = 24);

// Координаты станций
list<pair<double, double>> getStationList(vector<vector<Ray>> data);


void solveSle(Grid &grid, const vector<VectorSparse> &matrix, const vector<float> integrals, double error, bool onlyPositive = true);
void computeParametrs(Grid &crude, Grid &accurate, vector<VectorSparse> sleMatrix, vector<float> integrals, bool useSecondGrid, ElectronDensityDistribution *model, Dimension latitude, Dimension longitude, Dimension time, int intervals, int intervalsTime, double initialResidual);

double degreeToRadian(double degree);
double radianToDegree(double radian);

list<int> createListOfIntervals(int first, int last);

#endif
