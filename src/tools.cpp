#include <fstream>
#include <experimental/filesystem>
#include <fmt/format.h>
#include <fmt/ostream.h>

#include "tools.h"
#include "math.h"
#include "integration.h"

namespace fs = std::experimental::filesystem;

extern std::string pathToProcessedData;

void data_to_sle(const std::vector<std::vector<Ray>> &data, std::vector<VectorSparse> &phi, std::vector<double> &integrals, const Grid &test) {
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

void data_to_sle(const std::vector<std::vector<Ray>> &data, std::vector<VectorSparse> &phi, const Grid &test) {
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

double compute_residual(const Grid &x, const std::vector<VectorSparse> &A, const std::vector<double> &m) {
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

std::vector<double> compute_vector_residual(const Grid &x, const std::vector<VectorSparse> &A, const std::vector<double> &m) {
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

std::vector<std::vector<Ray>> get_data(const std::string &path, const unsigned startTime, const unsigned finishTime) {
    std::vector<std::vector<Ray>> data;
    if (fs::exists(path) && fs::is_directory(path)) {
        for (const auto& file : fs::directory_iterator(path)) {
            if (file.path().extension() == ".dat") {
                std::ifstream gps(file.path().string());
                if (!gps) {
                    fmt::print("Can't open file {}\n", file.path().string());
                } else {
                    std::vector<Ray> bundle;
                    int numberOfRays;
                    gps >> numberOfRays;

                    while (numberOfRays--) {
                        Ray ray;
                        gps >> ray;
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

std::list<std::pair<double, double>> get_station_list(std::vector<std::vector<Ray>> data) {
    std::list<std::pair<double, double>> stations;
    ChepmanLayer chepmanLayer;
    chepmanLayer.coordinateTransformation = std::make_unique<DecartToGeographic>();
    for (const auto& i : data) {
        for (const auto& j : i) {
            point station = j.station;
            chepmanLayer.coordinateTransformation->forward(station);
            bool found = false;
            for (const auto& k : stations) {
                if (k.first == station.R[0] && k.second == station.R[1]) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                stations.push_back(std::make_pair(station.R[0], station.R[1]));
            }
        }
    }
    return stations;
}


void solve_sle(Grid &grid, const std::vector<VectorSparse> &matrix, const std::vector<double> &integrals, const double error, const bool onlyPositive) {
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

void compute_parametrs(Grid &crude, Grid &accurate, const std::vector<VectorSparse> &sleMatrix, const std::vector<double> &integrals, const bool useSecondGrid,
 ElectronDensityDistribution &model, Dimension latitude, Dimension longitude, Dimension time, unsigned intervals, unsigned intervalsTime, double initialResidual) {
    latitude.to_radian();
    longitude.to_radian();
    double reconstructionSum = 0;
    double modelSum = 0;
    unsigned density = 100; // Количество точек по оси, по которым строится область
    const unsigned timeStart = static_cast<unsigned>(time.left / 3600);
    const unsigned timeFinish = static_cast<unsigned>(time.right / 3600);
    for (unsigned t = timeStart; t < timeFinish; ++t) {
        for (unsigned i = 0; i <= density; ++i) {
            for (unsigned j = 0; j <= density; ++j) {
                const double phi = latitude.left + (latitude.right - latitude.left) / density * i;
                const double theta = longitude.left + (longitude.right - longitude.left) / density * j;
                const double ray_time = t * 3600;

                point station(phi, theta, Re);
                point satellite(phi, theta, Re + 1000);
                model.coordinateTransformation->backward(station);
                model.coordinateTransformation->backward(satellite);

                Ray L(station, satellite, ray_time);

                Rectangle integral;

                const double crudeValue = crude(phi, theta, ray_time);
                const double modelValue = integral(L, model);

                modelSum += modelValue * modelValue;
                if (useSecondGrid) {
                    const double sumValue = crudeValue + accurate(phi, theta, ray_time);
                    reconstructionSum += (sumValue - modelValue) * (sumValue - modelValue);
                } else {
                    reconstructionSum += (crudeValue - modelValue) * (crudeValue - modelValue);
                }
            }
        }
    }

    latitude.to_degrees();
    longitude.to_degrees();

    std::ofstream parametrs(pathToProcessedData + "parametrs.txt", std::ios::app);
    auto residual = compute_residual(useSecondGrid ? accurate : crude, sleMatrix, integrals);
    fmt::print(parametrs, "{}\t{}\t{}\t{}\t{}\t{}\n",
               intervals,
               intervalsTime,
               residual / initialResidual,
               sqrt(reconstructionSum / modelSum),
               longitude.length() / intervals,
               time.length() / intervalsTime / 60);
    parametrs.close();
}

std::list<unsigned> create_intervals(unsigned first, unsigned last) {
    std::list<unsigned> foo;
    do {
        foo.push_back(first);
        first *= 2;
    } while (first <= last);
    return foo;
}
