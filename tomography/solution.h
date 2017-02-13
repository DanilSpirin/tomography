#pragma once
#include "grid.h"
#include "tools.h"

class Solution {
    ElectronDensityDistribution *model;
    vector<Grid> grids;
    vector<vector<Ray>> data;
    float latitudeLeft, longitudeLeft, timeLeft;
    float latitudeRight, longitudeRight, timeRight;
public:
    Solution();
    void setLimits(float latitudeLeft, float latitudeRight, float longitudeLeft, float longitudeRight, float timeLeft, float timeRight);
    void setModel(ElectronDensityDistribution model);
    void addGrid(int spaceIntervals, int timeIntervals);
    void addData(vector<vector<Ray>> data);
    void find();
    void print();
};