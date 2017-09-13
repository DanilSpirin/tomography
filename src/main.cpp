#include <fstream>
#include <cmath>
#include <fmt/format.h>

#include "reconstruction.h"
#include "integration.h"
#include "tools.h"
#include "limits.h"

std::string pathToData = "/home/ds/science/2003_302_1/";
std::string pathToProcessedData = pathToData + "tec_processed_model/";

int main(int argc, const char * argv[]) {

    const bool useSecondGrid = true;
    const bool calcParametrs = false;
    
    ChepmanLayer chepmanLayer;

    chepmanLayer.add_spot(point(30, 55), 12, pi / 24, 150, 0.5);
    chepmanLayer.add_spot(point(2.5, 45), 12, pi / 6, 300, 0.2);
    chepmanLayer.add_spot(point(12, 40), 12, pi / 4, 300, 0.2);

    chepmanLayer.add_wave(point(10, 40, Re), 32400, 1500, 0.3); // 25 minutes
    chepmanLayer.add_wave(point(20, 60, Re), 54000, 2100, 0.3); // 35 minutes
    chepmanLayer.add_wave(point(10, 40, Re), 32400, 1800, 0.3); // 20 minutes

    const unsigned timeStart = 0;
    const unsigned timeFinish = 24;

    std::list<unsigned> crudeIntervalsDim = createListOfIntervals(8, 8);
    std::list<unsigned> crudeIntervalsTime = createListOfIntervals(timeFinish - timeStart, timeFinish - timeStart);
    std::list<unsigned> accIntervalsDim = createListOfIntervals(36, 36);
    std::list<unsigned> accIntervalsTime = createListOfIntervals(240, 240);

    chepmanLayer.coordinateTransformation = std::make_unique<DecartToGeographic>();
    auto data = get_data(pathToData, timeStart, timeFinish);
    const auto stations = getStationList(data);

    std::ofstream station_file(pathToProcessedData + "stations_check.txt");
    for (const auto &[lat, lon] : stations) {
        station_file << radianToDegree(lat) << ' ' << radianToDegree(lon) << '\n';
    }
    station_file.close();

    Rectangle integral;
    for (auto &i : data) {
        for (auto &j : i) {
            j.integral = integral(j, chepmanLayer);
        }
    }

    Dimension latitude, longitude, time;

    Grid crude, accurate;

    if (calcParametrs) {
        std::ofstream parametrs(pathToProcessedData + "parametrs.txt");
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
                            computeParametrs(crude, accurate, accurateSleMatrix, accurateIntegrals, true,
                                chepmanLayer, latitude, longitude, time, intervalsAcc, intervalsTime, initialResidual);
                        }
                    }
                }
            } else if (calcParametrs) {
                computeParametrs(crude, accurate, crudeSleMatrix, crudeIntegrals, false,
                    chepmanLayer, latitude, longitude, time, intervals, intervalsT, initialResidual);
            }
        }

    }

    if (!calcParametrs) {

        std::ofstream crudeOut, accurateOut, sumOut, modelOut;

        Limits crude_limits, accurate_limits, model_limits, sum_limits;

        latitude.toRadian();
        longitude.toRadian();

        for (unsigned t = timeStart; t < timeFinish + 1; ++t) {
            const int density = 100;
            auto path = fmt::format("{}{}{:02}{}", pathToProcessedData, "time_", t, "first.txt");
            crudeOut.open(path);
            path = fmt::format("{}{}{:02}{}", pathToProcessedData, "time_", t, "model.txt");
            modelOut.open(path);
            if (useSecondGrid) {
                path = fmt::format("{}{}{:02}{}", pathToProcessedData, "time_", t, "second.txt");
                accurateOut.open(path);
                path = fmt::format("{}{}{:02}{}", pathToProcessedData, "time_", t, "sum.txt");
                sumOut.open(path);
            }
            for (int i = 0; i <= density; ++i) {
                for (int j = 0; j <= density; ++j) {
                    const double phi = latitude.left + (latitude.right - latitude.left) / density * i;
                    const double theta = longitude.left + (longitude.right - longitude.left) / density * j;
                    double time = t * 3600;

                    point station(phi, theta, Re);
                    point satellite(phi, theta, Re + 1000);
                    chepmanLayer.coordinateTransformation->backward(station);
                    chepmanLayer.coordinateTransformation->backward(satellite);
                    Ray L(station, satellite, time);

                    const double crudeValue = crude(phi, theta, time);
                    const double modelValue = integral(L, chepmanLayer);

                    crude_limits.update(crudeValue);
                    model_limits.update(modelValue);

                    crudeOut << crudeValue;
                    (j != density) ? (crudeOut << " ") : (crudeOut << '\n');

                    modelOut << modelValue;
                    (j != density) ? (modelOut << " ") : (modelOut << '\n');

                    if (useSecondGrid) {
                        const double accurateValue = accurate(phi, theta, time);
                        const double sumValue = crudeValue + accurateValue;

                        accurate_limits.update(accurateValue);
                        sum_limits.update(sumValue);

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

        std::ofstream gridLimits(pathToProcessedData + "limits.txt");

        gridLimits << latitude << '\n' << longitude << '\n'
            << crude_limits << '\n' << model_limits << '\n';
        if (useSecondGrid) {
            gridLimits << accurate_limits << '\n' << sum_limits << '\n';
        }
        gridLimits.close();
    }

    return 0;
}
