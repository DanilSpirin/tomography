#include <iostream>
#include <fstream>

#include <dirent.h>

#include "tools.h"
#include "math.h"
#include "integration.h"
#include <cstring>

extern unsigned timeStart, timeFinish;
extern std::string pathToProcessedData;

void dataToSle(const std::vector<std::vector<Ray>> &data, std::vector<VectorSparse> &phi, std::vector<double> &integrals, const Grid &test) {
    phi.clear();
    integrals.clear();
    for (int j = 0; j < data.size(); ++j) {
        for (int i = 0; i < data[j].size() - 1; ++i) {
//          Разностная схема
            phi.push_back((1 / cos(data[j][i + 1].angle)) * test.basis(data[j][i+1].phi, data[j][i + 1].thetta, data[j][i + 1].time) - (1 / cos(data[j][i].angle)) * test.basis(data[j][i].phi, data[j][i].thetta, data[j][i].time));
            integrals.push_back(data[j][i + 1].integral - data[j][i].integral);
        }
    }
}

void dataToSle(const std::vector<std::vector<Ray>> &data, std::vector<VectorSparse> &phi, const Grid &test) {
    phi.clear();
    for (int j = 0; j < data.size(); ++j) {
        for (int i = 0; i < data[j].size() - 1; ++i) {
            phi.push_back((1 / cos(data[j][i + 1].angle)) * test.basis(data[j][i + 1].phi, data[j][i + 1].thetta, data[j][i + 1].time) - (1 / cos(data[j][i].angle)) * test.basis(data[j][i].phi, data[j][i].thetta, data[j][i].time));
        }
    }
}

double computeResidual(const Grid &x, const std::vector<VectorSparse> &A, const std::vector<double> &m) {
    std::vector<double> difference(m.size(), 0);
    for (int i = 0; i < A.size(); ++i) {
        for (int j = 0; j < A[i].getSize(); ++j) {
            difference[i] += A[i].getPhi(j) * x[A[i].getNumber(j)];
        }
        difference[i] -= m[i];
    }
    double sum = 0;
    for (int i = 0; i < difference.size(); ++i) {
        sum += difference[i] * difference[i];
    }

    return sqrt(sum);
}

std::vector<double> computeVectorResidual(const Grid &x, const std::vector<VectorSparse> &A, const std::vector<double> &m) {
    std::vector<double> difference(m.size(), 0);
    for (int i = 0; i < A.size(); ++i) {
        for (int j = 0; j < A[i].getSize(); ++j) {
            difference[i] += A[i].getPhi(j) * x[A[i].getNumber(j)];
        }
        difference[i] -= m[i];
    }
    for (int i = 0; i < difference.size(); ++i) {
        difference[i] = -difference[i];
    }
    return difference;
}

// отсеиваем скрытые файлы и файлы не являющиеся данными
int sel(const struct dirent *d) {
    if (!strcmp(".", d->d_name) || !strcmp("..", d->d_name) || !strcmp(".DS_Store", d->d_name)) {
        return 0;
    }
    char* p = strstr(d->d_name, ".dat");
    return p ? 1 : 0;
}
std::vector<std::vector<Ray>> getData(const char *pathToData, int startTime, int finishTime) {
    startTime *= 3600;
    finishTime *= 3600;
    dirent **nameList;
    int numberOfBundles = scandir(pathToData, &nameList, sel, 0);
    char temp[100];
    std::vector<std::vector<Ray>> data;

    if (numberOfBundles < 0) {
        perror("Scandir");
    }
    else if (numberOfBundles == 0) {
        std::cout << "No data files found!" << std::endl;
    }
    else {
        while (numberOfBundles--) {
            std::vector<Ray> tempBundle;
            sprintf(temp, "%s%s", pathToData, nameList[numberOfBundles]->d_name);
            std::ifstream gps(temp);

            if (!gps) {
                std::cout << "Can't open file " << temp << std::endl;
            }
            else {
                int numberOfRays;
                gps >> numberOfRays;

                while (numberOfRays--) {
                    Ray tempRay;
                    gps >> tempRay;
                    if (tempRay.time >= startTime && tempRay.time <= finishTime) {
                        tempRay.computeParameters();
                        tempBundle.push_back(tempRay);
                    }
                }
                gps.close();
            }
            if (tempBundle.size() > 1) {
                data.push_back(tempBundle);
            }
            free(nameList[numberOfBundles]);
        }
    }

    free(nameList);

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
    return stations;
}


void solveSle(Grid &grid, const std::vector<VectorSparse> &matrix, const std::vector<double> integrals, double error, bool onlyPositive) {
    double initialResidual = computeResidual(grid, matrix, integrals);
    double firstRes, secondRes, currentRes;
    double limit = 0;
    double iterations = 50;
    double counter = 0;
    for (int i = 0; i < iterations; ++i) {
        iterationSirt(grid, matrix, integrals, onlyPositive);
    }
    firstRes = computeResidual(grid, matrix, integrals) / initialResidual;
    for (int i = 0; i < iterations; ++i) {
        iterationSirt(grid, matrix, integrals, onlyPositive);
    }
    secondRes = computeResidual(grid, matrix, integrals) / initialResidual;
    limit = (iterations * 2 * secondRes - iterations * firstRes) / iterations;
    std::cout << limit << std::endl;
    currentRes = secondRes;
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


double degreeToRadian(double degree){
    return degree / 180 * pi;
}

double radianToDegree(double radian) {
    return radian * 180 / pi;
}

void computeParametrs(Grid &crude, Grid &accurate, std::vector<VectorSparse> sleMatrix, std::vector<double> integrals, bool useSecondGrid, ElectronDensityDistribution &model, Dimension latitude, Dimension longitude, Dimension time, int intervals, int intervalsTime, double initialResidual) {
    latitude.toRadian();
    longitude.toRadian();
    double reconstructionSum = 0, modelSum = 0;
    int density = 100; // Количество точек по оси, по которым строится область
    for (int t = timeStart; t < timeFinish; ++t) {
        for (int i = 0; i <= density; ++i) {
            for (int j = 0; j <= density; ++j) {
                double phi = latitude.left + (latitude.right - latitude.left) / density * i;
                double theta = longitude.left + (longitude.right - longitude.left) / density * j;
                double time = t * 3600;

                point station(phi, theta, Re);
                point satellite(phi, theta, Re + 1000);
                model.coordinateTransformation->backward(station);
                model.coordinateTransformation->backward(satellite);

                Ray L(station, satellite, time);

                double crudeValue, accurateValue, sumValue, modelValue;
                Rectangle integral;

                crudeValue = crude(phi, theta, time);
                modelValue = integral(L, model);

                modelSum += modelValue * modelValue;
                if (useSecondGrid) {
                    accurateValue = accurate(phi, theta, time);
                    sumValue = crudeValue + accurateValue;
                }
                if (useSecondGrid) {
                    reconstructionSum += (sumValue - modelValue) * (sumValue - modelValue);
                }
                else {
                    reconstructionSum += (crudeValue - modelValue) * (crudeValue - modelValue);
                }
            }
        }
    }

    latitude.toDegrees();
    longitude.toDegrees();

    std::ofstream parametrs((pathToProcessedData+"parametrs.txt").c_str(), std::ios::app);

    parametrs << intervals << '\t' << intervalsTime << '\t';
    if (useSecondGrid) {
        parametrs << computeResidual(accurate, sleMatrix, integrals) / initialResidual << '\t';
    }
    else {
        parametrs << computeResidual(crude, sleMatrix, integrals) / initialResidual << '\t';
    }
    parametrs << sqrt(reconstructionSum / modelSum) << '\t' << (longitude.right - longitude.left) / intervals << '\t' << (time.right - time.left) / intervalsTime / 60 << std::endl;
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
