#pragma once

#include "tools.h"

class Solution {
    ElectronDensityDistribution *model;
    std::vector<Grid> grids;
    std::vector<std::vector<Ray>> data;
    double latitudeLeft, longitudeLeft, timeLeft;
    double latitudeRight, longitudeRight, timeRight;
public:
    Solution() = default;
    ~Solution() = default;
    void set_limits(double latitudeLeft, double latitudeRight,
                    double longitudeLeft, double longitudeRight,
                    double timeLeft, double timeRight);
    void set_model(ElectronDensityDistribution &model);
    void add_grid(unsigned spaceIntervals, unsigned timeIntervals);
    void add_data(std::vector<std::vector<Ray>> &&data);
    void find();
    void print(const std::string& output_path);
};
