#include <fstream>

#include <fmt/format.h>
#include <fmt/ostream.h>

#include "limits.h"
#include "solution.h"


void Solution::set_limits(float latitudeLeft, float latitudeRight,
                          float longitudeLeft, float longitudeRight,
                          float timeLeft, float timeRight) {
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
    Dimension latitude(latitudeLeft, latitudeRight, spaceIntervals);
    Dimension longitude(longitudeLeft, longitudeRight, spaceIntervals);
    Dimension time(timeLeft * 3600, timeRight * 3600, timeIntervals);
    Grid grid(latitude, longitude, time);
    grids.push_back(grid);
}

void Solution::add_data(std::vector<std::vector<Ray>> &&_data) {
    this->data = _data;
}

void Solution::find(const Solver& solver) {
    auto numberOfGrids = grids.size();
    if (numberOfGrids > 0) {
        std::vector<float> currentIntegrals;
        std::vector<VectorSparse> currentSleMatrix;

        data_to_sle(data, currentSleMatrix, currentIntegrals, grids.at(0));
        solve_sle(grids.at(0), currentSleMatrix, currentIntegrals, 0.1f, solver);

        for (unsigned i = 1; i < numberOfGrids; ++i) {
            currentIntegrals = compute_vector_residual(grids.at(i - 1), currentSleMatrix, currentIntegrals);
            data_to_sle(data, currentSleMatrix, grids.at(i));
            solve_sle(grids.at(i), currentSleMatrix, currentIntegrals, 0.15f, solver, false);
        }

    } else {
        fmt::print("No grids, can't solve\n");
    }
}

void Solution::print(const std::string& output_path) {
    Limits limit;
    Dimension latitude(latitudeLeft, latitudeRight, 0);
    Dimension longitude(longitudeLeft, longitudeRight, 0);

    latitude.to_radian();
    longitude.to_radian();

    const unsigned density = 100;

    for (unsigned i = static_cast<unsigned>(timeLeft); i < timeRight + 1; ++i) {
        auto path = fmt::format("{}{}{:02}{}", output_path, "time_", i, ".txt");
        std::ofstream out(path);
        for (unsigned x = 0; x <= density; ++x) {
            for (unsigned y = 0; y <= density; ++y) {
                const float phi = latitude.left + (latitude.right - latitude.left) / density * x;
                const float theta = longitude.left + (longitude.right - longitude.left) / density * y;
                const float time = i * 3600;

                float sum = 0;
                for (const auto &item : grids) {
                    sum += item(phi, theta, time);
                }

                limit.update(sum);
                fmt::print(out, "{:.2f}{}", sum, (y != density) ? ' ' : '\n');
            }
        }
        out.close();
    }

    latitude.to_degrees();
    longitude.to_degrees();

    std::ofstream limits(output_path + "limits.txt");
    fmt::print(limits, "{}\n{}\n{}", latitude, longitude, limit);
    limits.close();
}
