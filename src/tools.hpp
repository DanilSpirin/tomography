#pragma once

#include <set>
#include <vector>

#include "distribution.hpp"
#include "grid.hpp"
#include "ray.hpp"
#include "reconstruction.hpp"
#include "vector_sparse.hpp"

using SparseMatrix = std::vector<VectorSparse>;
using SleMatrix = std::vector<std::vector<Ray>>;

// Forming task SLE
void data_to_sle(const SleMatrix &data, SparseMatrix &phi, std::vector<float> &integrals, const Grid &test);
void data_to_sle(const SleMatrix &data, SparseMatrix &phi, const Grid &test);

// Calculating residual
float compute_residual(const Grid &x, const SparseMatrix &A, const std::vector<float> &m);
std::vector<float> compute_vector_residual(const Grid &x, const SparseMatrix &A, const std::vector<float> &m);

// Reading data from specified path
SleMatrix get_data(const std::string &path, const unsigned startTime = 0, const unsigned finishTime = 24);

// Station coordinates
std::set<std::pair<float, float>> get_stations(const SleMatrix& data);

void solve_sle(Grid &grid, const SparseMatrix &matrix, const std::vector<float> &integrals, const float error, const Solver &solver, const bool onlyPositive = true);

float degree_to_radian(const float degree);
float radian_to_degree(const float radian);

std::vector<unsigned> create_intervals(const unsigned first, const unsigned last);
