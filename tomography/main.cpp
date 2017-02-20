#include <iostream>
#include <fstream>
#include <cmath>
#include <climits>

#include "reconstruction.h"
#include "integration.h"
#include "tools.h"

bool useSecondGrid = true;
bool calcParametrs = false;

unsigned timeStart, timeFinish;

std::string pathToData = "/Users/imaginary/Documents/Science/2003_302_1/";
std::string pathToProcessedData = pathToData + "tec_processed_model/";

int main(int argc, const char * argv[]) {
    std::ofstream crudeOut, accurateOut, sumOut, modelOut;

    Rectangle integral;
    ChepmanLayer chepmanLayer;

    chepmanLayer.add_spot(point(30, 55), 12, pi / 24, 150, 0.5);
    chepmanLayer.add_spot(point(2.5, 45), 12, pi / 6, 300, 0.2);
    chepmanLayer.add_spot(point(12, 40), 12, pi / 4, 300, 0.2);

    chepmanLayer.add_wave(point(10, 40, Re), 32400, 1500, 0.3); // 25 minutes
    chepmanLayer.add_wave(point(20, 60, Re), 54000, 2100, 0.3); // 35 minutes
    chepmanLayer.add_wave(point(10, 40, Re), 32400, 1800, 0.3); // 20 minutes

    timeStart = 0;
    timeFinish = 24;

    std::list<unsigned> crudeIntervalsDim = createListOfIntervals(8, 8);
    std::list<unsigned> crudeIntervalsTime = createListOfIntervals(timeFinish - timeStart, timeFinish - timeStart);
    std::list<unsigned> accIntervalsDim = createListOfIntervals(36, 36);
    std::list<unsigned> accIntervalsTime = createListOfIntervals(240, 240);

    chepmanLayer.coordinateTransformation = new DecartToGeographic;

    std::vector<std::vector<Ray>> data = getData(pathToData.c_str(), timeStart, timeFinish);

    std::list<std::pair<double, double>> stations = getStationList(data);

    std::ofstream station_file((pathToProcessedData + "stations_check.txt").c_str());
    for (const auto &st : stations) {
        station_file << radianToDegree(st.first) << ' ' << radianToDegree(st.second) << '\n';
    }
    station_file.close();

    for (auto &i : data) {
        for (auto &j : i) {
            j.integral = integral(j, chepmanLayer);
        }
    }

    Dimension latitude, longitude, time;

    Grid crude, accurate;

    if (calcParametrs) {
        std::ofstream parametrs((pathToProcessedData + "parametrs.txt").c_str());
        parametrs.close();
    }

    for (unsigned intervals : crudeIntervalsDim) {
        for (unsigned intervalsT : crudeIntervalsTime) {
            latitude = Dimension(-10.0, 40.0, intervals);
            longitude = Dimension(30.0, 70.0, intervals);
            time = Dimension(double(timeStart * 3600), double(timeFinish * 3600), intervalsT);

            std::vector<VectorSparse> crudeSleMatrix;
            std::vector<double> crudeIntegrals;

            crude.set(latitude, longitude, time);
            dataToSle(data, crudeSleMatrix, crudeIntegrals, crude);
            double initialResidual = computeResidual(crude, crudeSleMatrix, crudeIntegrals);

            solveSle(crude, crudeSleMatrix, crudeIntegrals, 0.15);

            if (useSecondGrid) {
                std::vector<double> accurateIntegrals;
                accurateIntegrals = computeVectorResidual(crude, crudeSleMatrix, crudeIntegrals);
                std::vector<VectorSparse>().swap(crudeSleMatrix);
                std::vector<double>().swap(crudeIntegrals);
                for (unsigned intervalsAcc : accIntervalsDim) {
                    for (unsigned intervalsTime : accIntervalsTime) {
                        latitude = Dimension(-10.0, 40.0, intervalsAcc);
                        longitude = Dimension(30.0, 70.0, intervalsAcc);
                        time = Dimension(double(timeStart * 3600), double(timeFinish * 3600), intervalsTime);

                        std::vector<VectorSparse> accurateSleMatrix;

                        accurate.set(latitude, longitude, time);
                        dataToSle(data, accurateSleMatrix, accurate);

                        initialResidual = computeResidual(accurate, accurateSleMatrix, accurateIntegrals);

                        solveSle(accurate, accurateSleMatrix, accurateIntegrals, 0.1, false);

                        if (calcParametrs) {
                            computeParametrs(crude, accurate, accurateSleMatrix, accurateIntegrals, true, chepmanLayer, latitude, longitude, time, intervalsAcc, intervalsTime, initialResidual);
                        }
                    }
                }
            }
            else if (calcParametrs) {
                computeParametrs(crude, accurate, crudeSleMatrix, crudeIntegrals, false, chepmanLayer, latitude, longitude, time, intervals, intervalsT, initialResidual);
            }
        }

    }

    if (!calcParametrs) {

        double crudeMin = INT_MAX, crudeMax = INT_MIN;
        double accurateMin = INT_MAX, accurateMax = INT_MIN;
        double sumMin = INT_MAX, sumMax = INT_MIN;
        double modelMin = INT_MAX, modelMax = INT_MIN;

        latitude.toRadian();
        longitude.toRadian();

        char path[100];
        int density = 100;
        for (unsigned t = timeStart; t < timeFinish + 1; ++t) {
            sprintf(path, "%s%s%02d%s", pathToProcessedData.c_str(), "time_", t, "first.txt");
            crudeOut.open(path);
            sprintf(path, "%s%s%02d%s", pathToProcessedData.c_str(), "time_", t, "model.txt");
            modelOut.open(path);
            if (useSecondGrid) {
                sprintf(path, "%s%s%02d%s", pathToProcessedData.c_str(), "time_", t, "second.txt");
                accurateOut.open(path);
                sprintf(path, "%s%s%02d%s", pathToProcessedData.c_str(), "time_", t, "sum.txt");
                sumOut.open(path);
            }
            for (int i = 0; i <= density; ++i) {
                for (int j = 0; j <= density; ++j) {
                    double phi = latitude.left + (latitude.right - latitude.left) / density * i;
                    double theta = longitude.left + (longitude.right - longitude.left) / density * j;
                    double time = t * 3600;

                    point station(phi, theta, Re);
                    point satellite(phi, theta, Re + 1000);
                    chepmanLayer.coordinateTransformation->backward(station);
                    chepmanLayer.coordinateTransformation->backward(satellite);

                    Ray L(station, satellite, time);

                    double crudeValue, modelValue;

                    crudeValue = crude(phi, theta, time);
                    modelValue = integral(L, chepmanLayer);

                    if (crudeValue < crudeMin) {
                        crudeMin = crudeValue;
                    }
                    if (crudeValue > crudeMax) {
                        crudeMax = crudeValue;
                    }

                    if (modelValue < modelMin) {
                        modelMin = modelValue;
                    }
                    if (modelValue > modelMax) {
                        modelMax = modelValue;
                    }

                    crudeOut << crudeValue;
                    (j != density) ? (crudeOut << " ") : (crudeOut << '\n');

                    modelOut << modelValue;
                    (j != density) ? (modelOut << " ") : (modelOut << '\n');

                    if (useSecondGrid) {
                        double accurateValue, sumValue;
                        accurateValue = accurate(phi, theta, time);
                        sumValue = crudeValue + accurateValue;

                        if (accurateValue < accurateMin) {
                            accurateMin = accurateValue;
                        }
                        if (accurateValue > accurateMax) {
                            accurateMax = accurateValue;
                        }

                        if (sumValue < sumMin) {
                            sumMin = sumValue;
                        }
                        if (sumValue > sumMax) {
                            sumMax = sumValue;
                        }

                        accurateOut << accurateValue;
                        (j != density) ? (accurateOut << " ") : (accurateOut << '\n');
                        sumOut << sumValue;
                        (j != density) ? (sumOut << " ") : (sumOut << '\n');
                    }
                }
            }
            crudeOut.close();
            modelOut.close();
            if (useSecondGrid) {
                accurateOut.close();
                sumOut.close();
            }
        }

        latitude.toDegrees();
        longitude.toDegrees();

        std::ofstream gridLimits((pathToProcessedData + "limits.txt").c_str());

        gridLimits << latitude << '\n' << longitude << '\n'
                   << floor(crudeMin) << ' ' << ceil(crudeMax) << '\n'
                   << floor(modelMin) << ' ' << ceil(modelMax) << '\n';
        if (useSecondGrid) {
            gridLimits  << floor(accurateMin) << ' ' << ceil(accurateMax) << '\n'
                        << floor(sumMin) << ' ' << ceil(sumMax) << '\n';
        }
        gridLimits.close();
    }

    delete chepmanLayer.coordinateTransformation;
    return 0;
}