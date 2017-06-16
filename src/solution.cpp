#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

#include "solution.h"
#include <climits>
extern std::string pathToProcessedData;

Solution::Solution() : grids(0){};

void Solution::setLimits(float latitudeLeft, float latitudeRight, float longitudeLeft, float longitudeRight, float timeLeft, float timeRight) {
    this->latitudeLeft = latitudeLeft;
    this->latitudeRight = latitudeRight;
    this->longitudeLeft = longitudeLeft;
    this->longitudeRight = longitudeRight;
    this->timeLeft = timeLeft;
    this->timeRight = timeRight;
}

void Solution::setModel(ElectronDensityDistribution &model) {
    this->model = &model;
}

void Solution::addGrid(int spaceIntervals, int timeIntervals) {
    Grid foo;
    Dimension latitude(latitudeLeft, latitudeRight, spaceIntervals * 2);
    Dimension longitude(longitudeLeft, longitudeRight, spaceIntervals);
    Dimension time(timeLeft * 3600, timeRight * 3600, timeIntervals);
    foo.set(latitude, longitude, time);
    grids.push_back(foo);
}

void Solution::addData(std::vector<std::vector<Ray>> _data) {
    this->data = _data;
}

void Solution::find() {
    int numberOfGrids = (int)grids.size();
    if (numberOfGrids > 0) {
        std::vector<double> currentIntegrals;
        std::vector<VectorSparse> currentSleMatrix;

        dataToSle(data, currentSleMatrix, currentIntegrals, grids.at(0));
        solveSle(grids.at(0), currentSleMatrix, currentIntegrals, 0.1);

        for (int i = 1; i < numberOfGrids; ++i) {
            currentIntegrals = computeVectorResidual(grids.at(i - 1), currentSleMatrix, currentIntegrals);
            dataToSle(data, currentSleMatrix, grids.at(i));
            solveSle(grids.at(i), currentSleMatrix, currentIntegrals, 0.15, false);
        }

    } else {
        std::cout << "No grids, can't solve" << std::endl;
    }
}

void Solution::print() {
    double min = INT_MAX, max = INT_MIN;
    Dimension latitude(latitudeLeft, latitudeRight, 0);
    Dimension longitude(longitudeLeft, longitudeRight, 0);

    latitude.toRadian();
    longitude.toRadian();

    int density = 150;

    std::ofstream out;
    char path[100];
    for (int i = (int)timeLeft; i < timeRight + 1; ++i) {
        sprintf(path, "%s%s%02d%s", pathToProcessedData.c_str(), "time_", i, ".txt");
        out.open(path);
        for (int x = 0; x <= density; ++x) {
            for (int y = 0; y <= density; ++y) {
                double phi = latitude.left + (latitude.right - latitude.left) / density * x;
                double theta = longitude.left + (longitude.right - longitude.left) / density * y;
                double time = i * 3600;

                double sum = 0;
                for (auto &&item : grids) {
                    sum += item(phi, theta, time);
                }
                if (sum < min) {
                    min = sum;
                }
                if (sum > max) {
                    max = sum;
                }

                out << sum;
                (y != density) ? (out << " ") : (out << std::endl);
            }
        }
        out.close();
    }

    latitude.toDegrees();
    longitude.toDegrees();

    sprintf(path, "%s%s", pathToProcessedData.c_str(), "limits.txt");
    std::ofstream limits(path);
    limits << latitude << std::endl << longitude << std::endl << floor(min) << ' ' << ceil(max);
    limits.close();
}