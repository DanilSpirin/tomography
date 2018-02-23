#pragma once

#include "tools.hpp"

class Solution {
    ElectronDensityDistribution *model;
    std::vector<Grid> grids;
    std::vector<std::vector<Ray>> data;
    float latitudeLeft, longitudeLeft, timeLeft;
    float latitudeRight, longitudeRight, timeRight;
public:
    Solution() = default;
    ~Solution() = default;
    void set_limits(float latitudeLeft, float latitudeRight,
                    float longitudeLeft, float longitudeRight,
                    float timeLeft, float timeRight);
    void set_model(ElectronDensityDistribution &model);
    void add_grid(unsigned spaceIntervals, unsigned timeIntervals);
    void add_data(std::vector<std::vector<Ray>> &&data);
    void find(const Solver &solver);
    void print(const std::string& output_path);
};
