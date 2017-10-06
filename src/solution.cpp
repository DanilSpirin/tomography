#include <fstream>
#include <cmath>
#include <string>
#include <numeric>
#include <fmt/format.h>

#include "solution.h"
#include "limits.h"

extern std::string pathToProcessedData;

void Solution::set_limits(double latitudeLeft, double latitudeRight, double longitudeLeft, double longitudeRight, double timeLeft, double timeRight) {
    this->latitudeLeft = latitudeLeft;
    this->latitudeRight = latitudeRight;
    this->longitudeLeft = longitudeLeft;
    this->longitudeRight = longitudeRight;
    this->timeLeft = timeLeft;
    this->timeRight = timeRight;
}

void Solution::set_model(ElectronDensityDistribution &model) {
    this->model = &model;
}

void Solution::add_grid(unsigned spaceIntervals, unsigned timeIntervals) {
    Grid foo;
    Dimension latitude(latitudeLeft, latitudeRight, spaceIntervals);
    Dimension longitude(longitudeLeft, longitudeRight, spaceIntervals);
    Dimension time(timeLeft * 3600, timeRight * 3600, timeIntervals);
    foo.set(latitude, longitude, time);
    grids.push_back(foo);
}

void Solution::add_data(std::vector<std::vector<Ray>> &&_data) {
    this->data = _data;
}

void Solution::find() {
    auto numberOfGrids = grids.size();
    if (numberOfGrids > 0) {
        std::vector<double> currentIntegrals;
        std::vector<VectorSparse> currentSleMatrix;

        data_to_sle(data, currentSleMatrix, currentIntegrals, grids.at(0));
        solve_sle(grids.at(0), currentSleMatrix, currentIntegrals, 0.1);

        for (unsigned i = 1; i < numberOfGrids; ++i) {
            currentIntegrals = compute_vector_residual(grids.at(i - 1), currentSleMatrix, currentIntegrals);
            data_to_sle(data, currentSleMatrix, grids.at(i));
            solve_sle(grids.at(i), currentSleMatrix, currentIntegrals, 0.15, false);
        }

    } else {
        fmt::print("No grids, can't solve\n");
    }
}

void Solution::print() {
    Limits limit;
    Dimension latitude(latitudeLeft, latitudeRight, 0);
    Dimension longitude(longitudeLeft, longitudeRight, 0);

    latitude.to_radian();
    longitude.to_radian();

    const unsigned density = 100;

    for (unsigned i = static_cast<unsigned>(timeLeft); i < timeRight + 1; ++i) {
        auto path = fmt::format("{}{}{:02}{}", pathToProcessedData, "time_", i, ".txt");
        std::ofstream out(path);
        for (unsigned x = 0; x <= density; ++x) {
            for (unsigned y = 0; y <= density; ++y) {
                const double phi = latitude.left + (latitude.right - latitude.left) / density * x;
                const double theta = longitude.left + (longitude.right - longitude.left) / density * y;
                const double time = i * 3600;

                double sum = 0;
                for (const auto &item : grids) {
                    sum += item(phi, theta, time);
                }

                limit.update(sum);

                out << sum;
                (y != density) ? (out << " ") : (out << std::endl);
            }
        }
        out.close();
    }

    latitude.to_degrees();
    longitude.to_degrees();

    std::ofstream limits(pathToProcessedData + "limits.txt");
    limits << latitude << '\n' << longitude << '\n' << limit;
    limits.close();
}
