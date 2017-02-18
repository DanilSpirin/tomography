#pragma once

#include <vector>
#include "vector_sparse.h"
#include "grid.h"

void iterationArt(Grid &grid, const std::vector<VectorSparse> &a, const std::vector<double> &m, bool onlyPositive=true);
void iterationSirt(Grid &grid, const std::vector<VectorSparse> &a, const std::vector<double> &m, bool onlyPositive=true);