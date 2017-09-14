#include <fstream>
#include <cmath>

#include "integration.h"
#include "tools.h"
#include "solution.h"

std::string pathToData = "/home/ds/science/2003_302_1/";
std::string pathToProcessedData = pathToData + "tec_processed_model/";

int main(int argc, const char * argv[]) {

    ChepmanLayer chepmanLayer;

    chepmanLayer.add_spot(point(30, 55), 12, pi / 24, 150, 0.5);
    chepmanLayer.add_spot(point(2.5, 45), 12, pi / 6, 300, 0.2);
    chepmanLayer.add_spot(point(12, 40), 12, pi / 4, 300, 0.2);

    chepmanLayer.add_wave(point(10, 40, Re), 32400, 1500, 0.3); // 25 minutes
    chepmanLayer.add_wave(point(20, 60, Re), 54000, 2100, 0.3); // 35 minutes
    chepmanLayer.add_wave(point(10, 40, Re), 32400, 1800, 0.3); // 20 minutes

    const unsigned timeStart = 0;
    const unsigned timeFinish = 24;

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

    Solution solution;
    solution.set_limits(-10.0, 40.0, 30.0, 70.0, timeStart, timeFinish);
    solution.add_data(std::move(data));
    solution.add_grid(8, 24);
    solution.add_grid(36, 240);
    solution.find();
    solution.print();

    return 0;
}
