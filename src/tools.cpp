#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>

#include "tools.h"
#include "math.h"
#include "integration.h"

extern std::string pathToProcessedData;

void dataToSle(const std::vector<std::vector<Ray>> &data, std::vector<VectorSparse> &phi, std::vector<double> &integrals, const Grid &test) {
    phi.clear();
    integrals.clear();
    for (int j = 0; j < data.size(); ++j) {
        for (int i = 0; i < data[j].size() - 1; ++i) {
            const auto left = test.basis(data[j][i + 1].phi, data[j][i + 1].thetta, data[j][i + 1].time) / cos(data[j][i + 1].angle);
            const auto right = test.basis(data[j][i].phi, data[j][i].thetta, data[j][i].time) / cos(data[j][i].angle);
            phi.push_back(left - right);
            integrals.push_back(data[j][i + 1].integral - data[j][i].integral);
        }
    }
}

void dataToSle(const std::vector<std::vector<Ray>> &data, std::vector<VectorSparse> &phi, const Grid &test) {
    phi.clear();
    for (int j = 0; j < data.size(); ++j) {
        for (int i = 0; i < data[j].size() - 1; ++i) {
            const auto left = test.basis(data[j][i + 1].phi, data[j][i + 1].thetta, data[j][i + 1].time) / cos(data[j][i + 1].angle);
            const auto right = test.basis(data[j][i].phi, data[j][i].thetta, data[j][i].time) / cos(data[j][i].angle);
            phi.push_back(left - right);
        }
    }
}

double computeResidual(const Grid &x, const std::vector<VectorSparse> &A, const std::vector<double> &m) {
    double sum = 0;
    for (int i = 0; i < A.size(); ++i) {
        double difference = 0;
        for (int j = 0; j < A[i].getSize(); ++j) {
            difference += A[i].getPhi(j) * x[A[i].getNumber(j)];
        }
        difference -= m[i];
        sum += (difference * difference);
    }
    return sqrt(sum);
}

std::vector<double> computeVectorResidual(const Grid &x, const std::vector<VectorSparse> &A, const std::vector<double> &m) {
    std::vector<double> difference(m.size(), 0);
    for (int i = 0; i < A.size(); ++i) {
        for (int j = 0; j < A[i].getSize(); ++j) {
            difference[i] -= A[i].getPhi(j) * x[A[i].getNumber(j)];
        }
        difference[i] += m[i];
    }
    return difference;
}

std::vector<std::vector<Ray>> get_data(const std::string &path, const unsigned startTime, const unsigned finishTime) {
    boost::filesystem::path p(path);
    std::vector<std::vector<Ray>> data;
    if (boost::filesystem::exists(p) && boost::filesystem::is_directory(p)) {
        for (const auto& file : boost::filesystem::directory_iterator(p)) {
            if (file.path().extension() == ".dat") {
                std::ifstream gps(file.path().string());
                if (!gps) {
                    std::cout << "Can't open file " << file.path().string() << std::endl;
                } else {
                    std::vector<Ray> bundle;
                    int numberOfRays;
                    gps >> numberOfRays;

                    while (numberOfRays--) {
                        Ray ray;
                        gps >> ray;
                        if (ray.time >= startTime * 3600 && ray.time <= finishTime * 3600) {
                            ray.computeParameters();
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

std::list<std::pair<double, double>> getStationList(std::vector<std::vector<Ray>> data) {
    std::list<std::pair<double, double>> stations;
    bool found;
    ChepmanLayer chepmanLayer;
    chepmanLayer.coordinateTransformation = new DecartToGeographic;
    for (const auto& i : data) {
        for (const auto& j : i) {
            point station = j.station;
            chepmanLayer.coordinateTransformation->forward(station);
            found = false;
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
    delete chepmanLayer.coordinateTransformation;
    return stations;
}


void solveSle(Grid &grid, const std::vector<VectorSparse> &matrix, const std::vector<double> &integrals, const double error, const bool onlyPositive) {
    const double initialResidual = computeResidual(grid, matrix, integrals);
    double iterations = 50;
    for (int i = 0; i < iterations; ++i) {
        iterationSirt(grid, matrix, integrals, onlyPositive);
    }
    const double firstRes = computeResidual(grid, matrix, integrals) / initialResidual;
    for (int i = 0; i < iterations; ++i) {
        iterationSirt(grid, matrix, integrals, onlyPositive);
    }
    const double secondRes = computeResidual(grid, matrix, integrals) / initialResidual;
    const double limit = (iterations * 2 * secondRes - iterations * firstRes) / iterations;
    std::cout << limit << std::endl;
    double currentRes = secondRes;
    unsigned counter = 0;
    while (currentRes / limit > 1 + error) {
        iterationSirt(grid, matrix, integrals, onlyPositive);
        currentRes = computeResidual(grid, matrix, integrals) / initialResidual;
        ++counter;
        if (counter > 500) {
            std::cout << "stopped at current/limit = " << currentRes / limit << std::endl;
            break;
        }
        std::cout << counter << std::endl;
    }
}


double degreeToRadian(const double degree) {
    return degree / 180 * pi;
}

double radianToDegree(const double radian) {
    return radian * 180 / pi;
}

void computeParametrs(Grid &crude, Grid &accurate, const std::vector<VectorSparse> &sleMatrix, const std::vector<double> &integrals, const bool useSecondGrid,
 ElectronDensityDistribution &model, Dimension latitude, Dimension longitude, Dimension time, int intervals, int intervalsTime, double initialResidual) {
    latitude.toRadian();
    longitude.toRadian();
    double reconstructionSum = 0;
    double modelSum = 0;
    int density = 100; // Количество точек по оси, по которым строится область
    const unsigned timeStart = time.left / 3600;
    const unsigned timeFinish = time.right / 3600;
    for (int t = timeStart; t < timeFinish; ++t) {
        for (int i = 0; i <= density; ++i) {
            for (int j = 0; j <= density; ++j) {
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

    latitude.toDegrees();
    longitude.toDegrees();

    std::ofstream parametrs(pathToProcessedData + "parametrs.txt", std::ios::app);

    parametrs << intervals << '\t' << intervalsTime << '\t';
    if (useSecondGrid) {
        parametrs << computeResidual(accurate, sleMatrix, integrals) / initialResidual << '\t';
    } else {
        parametrs << computeResidual(crude, sleMatrix, integrals) / initialResidual << '\t';
    }
    parametrs << sqrt(reconstructionSum / modelSum) << '\t' << longitude.length() / intervals << '\t' << time.length() / intervalsTime / 60 << std::endl;
    parametrs.close();

}

std::list<unsigned> createListOfIntervals(unsigned first, unsigned last) {
    std::list<unsigned> foo;
    do {
        foo.push_back(first);
        first *= 2;
    } while (first <= last);
    return foo;
}
