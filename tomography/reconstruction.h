#pragma once

#include <vector>
#include "vector_sparse.h"
#include "grid.h"
using namespace std;

void iterationArt(Grid &grid, const vector<VectorSparse> &a, const vector<double> &m, bool onlyPositive=true);
void iterationSirt(Grid &grid, const vector<VectorSparse> &a, const vector<double> &m, bool onlyPositive=true);