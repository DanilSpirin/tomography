#pragma once

#include <vector>

#include "grid.h"
#include "vector_sparse.h"

void iterationArt(Grid &grid, const std::vector<VectorSparse> &a, const std::vector<double> &m, bool onlyPositive=true);
void iterationSirt(Grid &grid, const std::vector<VectorSparse> &a, const std::vector<double> &m, bool onlyPositive=true);
