#include <fstream>
#include <experimental/filesystem>

#include <fmt/format.h>
#include <fmt/ostream.h>

#include "tools.h"

namespace fs = std::experimental::filesystem;

extern std::string pathToProcessedData;

using SparseMatrix = std::vector<VectorSparse>;
using SleMatrix = std::vector<std::vector<Ray>>;

void data_to_sle(const SleMatrix &data, SparseMatrix &phi, std::vector<double> &integrals, const Grid& test) {
    phi.clear();
    integrals.clear();
    for (const auto& it : data) {
        auto curr = std::cbegin(it);
        auto next = std::next(curr);
        for ( ; next != std::cend(it); ++curr, ++next) {
            const auto left  = test.basis(next->phi, next->thetta, next->time) / cos(next->angle);
            const auto right = test.basis(curr->phi, curr->thetta, curr->time) / cos(curr->angle);
            phi.push_back(left - right);
            integrals.push_back(next->integral - curr->integral);
        }
    }
}

void data_to_sle(const SleMatrix &data, SparseMatrix &phi, const Grid& test) {
    phi.clear();
    for (const auto& it : data) {
        auto curr = std::cbegin(it);
        auto next = std::next(curr);
        for ( ; next != std::cend(it); ++curr, ++next) {
            const auto left  = test.basis(next->phi, next->thetta, next->time) / cos(next->angle);
            const auto right = test.basis(curr->phi, curr->thetta, curr->time) / cos(curr->angle);
            phi.push_back(left - right);
        }
    }
}

double compute_residual(const Grid &x, const SparseMatrix &A, const std::vector<double> &m) {
    double sum = 0;
    for (unsigned i = 0; i < A.size(); ++i) {
        double difference = 0;
        for (unsigned j = 0; j < A[i].size(); ++j) {
            const auto element = A[i][j];
            difference += element.value * x[element.index];
        }
        difference -= m[i];
        sum += (difference * difference);
    }
    return sqrt(sum);
}

std::vector<double> compute_vector_residual(const Grid &x, const SparseMatrix &A, const std::vector<double> &m) {
    std::vector<double> difference(m.size(), 0);
    for (unsigned i = 0; i < A.size(); ++i) {
        for (unsigned j = 0; j < A[i].size(); ++j) {
            const auto element = A[i][j];
            difference[i] -= element.value * x[element.index];
        }
        difference[i] += m[i];
    }
    return difference;
}

SleMatrix get_data(const std::string &path, const unsigned startTime, const unsigned finishTime) {
    SleMatrix data;
    if (fs::exists(path) && fs::is_directory(path)) {
        for (const auto& file : fs::directory_iterator(path)) {
            if (file.path().extension() == ".dat") {
                std::ifstream gps(file.path().string());
                if (!gps) {
                    fmt::print("Can't open file {}\n", file.path().string());
                } else {
                    int numberOfRays;
                    gps >> numberOfRays;

                    std::vector<Ray> bundle;
                    bundle.reserve(numberOfRays);

                    Ray ray;
                    while (gps >> ray) {
                        if (ray.time >= startTime * 3600 && ray.time <= finishTime * 3600) {
                            ray.compute_parameters();
                            bundle.push_back(ray);
                        }
                    }
                    gps.close();
                    if (bundle.size() > 1) {
                        data.push_back(bundle);
                    }
                }
            }
        }
    }
    return data;
}

std::set<std::pair<double, double> > get_stations(const SleMatrix &data) {
    std::set<std::pair<double, double>> stations;
    ChepmanLayer chepmanLayer;
    chepmanLayer.coordinateTransformation = std::make_unique<DecartToGeographic>();
    for (const auto& i : data) {
        for (const auto& j : i) {
            point station = j.station;
            chepmanLayer.coordinateTransformation->forward(station);
            auto station_pair = std::make_pair(station.R[0], station.R[1]);
            if (stations.find(station_pair) == stations.end()) {
                stations.insert(station_pair);
            }
        }
    }
    return stations;
}

void solve_sle(Grid &grid, const SparseMatrix &matrix, const std::vector<double> &integrals,
        const double error, const bool onlyPositive) {
    const double initialResidual = compute_residual(grid, matrix, integrals);
    const double iterations = 50;
    for (int i = 0; i < iterations; ++i) {
        iterationSirt(grid, matrix, integrals, onlyPositive);
    }
    const double firstRes = compute_residual(grid, matrix, integrals) / initialResidual;
    for (int i = 0; i < iterations; ++i) {
        iterationSirt(grid, matrix, integrals, onlyPositive);
    }
    const double secondRes = compute_residual(grid, matrix, integrals) / initialResidual;
    const double limit = (iterations * 2 * secondRes - iterations * firstRes) / iterations;
    fmt::print("{}\n", limit);
    double currentRes = secondRes;
    unsigned counter = 0;
    while (currentRes / limit > 1 + error) {
        iterationSirt(grid, matrix, integrals, onlyPositive);
        currentRes = compute_residual(grid, matrix, integrals) / initialResidual;
        ++counter;
        if (counter > 500) {
            fmt::print("Stopped at current/limit = {}\n", currentRes / limit);
            break;
        }
        fmt::print("{}\n", counter);
    }
}

double degree_to_radian(const double degree) {
    return degree / 180 * pi;
}

double radian_to_degree(const double radian) {
    return radian * 180 / pi;
}

std::vector<unsigned> create_intervals(unsigned first, unsigned last) {
    std::vector<unsigned> foo;
    do {
        foo.push_back(first);
        first *= 2;
    } while (first <= last);
    return foo;
}
