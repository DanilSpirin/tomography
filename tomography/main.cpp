#include <iostream>
#include <fstream>
#include <cmath>

#include "reconstruction.h"
#include "integration.h"
#include "tools.h"

bool useSecondGrid = false;
bool calcParametrs = false;

int timeStart, timeFinish;

string pathToData = "/Users/imaginary/Documents/Science/2003_302_1/";
string pathToProcessedData = pathToData + "tec_processed_model/";

int main(int argc, const char * argv[]) {
    ofstream crudeOut, accurateOut, sumOut, modelOut;

    Rectangle integral;
    ChepmanLayer chepmanLayer;

    chepmanLayer.addSpot(point(30, 55), 12, pi / 24, 150, 0.5);
    chepmanLayer.addSpot(point(2.5, 45), 12, pi / 6, 300, 0.2);
    chepmanLayer.addSpot(point(12, 40), 12, pi / 4, 300, 0.2);

    chepmanLayer.addWave(point(10, 40, Re), 32400, 1500, 0.3); // 25 minutes
    chepmanLayer.addWave(point(20, 60, Re), 54000, 2100, 0.3); // 35 minutes
    chepmanLayer.addWave(point(10, 40, Re), 32400, 1800, 0.3); // 20 minutes

    timeStart = 0;
    timeFinish = 24;

    list<int> crudeIntervalsDim = createListOfIntervals(8, 8);
    list<int> crudeIntervalsTime = createListOfIntervals(timeFinish - timeStart, timeFinish - timeStart);
    list<int> accIntervalsDim = createListOfIntervals(36, 36);
    list<int> accIntervalsTime = createListOfIntervals(240, 240);

    chepmanLayer.coordinateTransformation = new DecartToGeographic;

    vector<vector<Ray>> data = getData(pathToData.c_str(), timeStart, timeFinish);

    list<pair<double, double>> stations = getStationList(data);

    ofstream station((pathToProcessedData + "stations_check.txt").c_str());
    for (const auto &st : stations) {
        station << radianToDegree(st.first) << ' ' << radianToDegree(st.second) << '\n';
    }
    station.close();

    for (auto &i : data) {
        for (auto &j : i) {
            j.integral = integral(j, chepmanLayer);
        }
    }

    Dimension latitude, longitude, time;

    Grid crude, accurate;

    if (calcParametrs) {
        ofstream parametrs((pathToProcessedData + "parametrs.txt").c_str());
        parametrs.close();
    }

    for (int intervals : crudeIntervalsDim) {
        for (int intervalsT : crudeIntervalsTime) {
            latitude = Dimension(-10.0, 40.0, intervals);
            longitude = Dimension(30.0, 70.0, intervals);
            time = Dimension(double(timeStart * 3600), double(timeFinish * 3600), intervalsT);

            vector<VectorSparse> crudeSleMatrix;
            vector<double> crudeIntegrals;

            crude.set(latitude, longitude, time);
            dataToSle(data, crudeSleMatrix, crudeIntegrals, crude);
            double initialResidual = computeResidual(crude, crudeSleMatrix, crudeIntegrals);

            solveSle(crude, crudeSleMatrix, crudeIntegrals, 0.15);

            if (useSecondGrid) {
                vector<double> accurateIntegrals;
                accurateIntegrals = computeVectorResidual(crude, crudeSleMatrix, crudeIntegrals);
                vector<VectorSparse>().swap(crudeSleMatrix);
                vector<double>().swap(crudeIntegrals);
                for (int intervalsAcc : accIntervalsDim) {
                    for (int intervalsTime : accIntervalsTime) {
                        latitude = Dimension(-10.0, 40.0, intervalsAcc);
                        longitude = Dimension(30.0, 70.0, intervalsAcc);
                        time = Dimension(double(timeStart * 3600), double(timeFinish * 3600), intervalsTime);

                        vector<VectorSparse> accurateSleMatrix;

                        accurate.set(latitude, longitude, time);
                        dataToSle(data, accurateSleMatrix, accurate);

                        initialResidual = computeResidual(accurate, accurateSleMatrix, accurateIntegrals);

                        solveSle(accurate, accurateSleMatrix, accurateIntegrals, 0.1, false);

                        if (calcParametrs) {
                            computeParametrs(crude, accurate, accurateSleMatrix, accurateIntegrals, true,chepmanLayer, latitude, longitude, time, intervalsAcc, intervalsTime, initialResidual);
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
        for (int t = timeStart; t < timeFinish + 1; ++t) {
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

        ofstream gridLimits((pathToProcessedData + "limits.txt").c_str());

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