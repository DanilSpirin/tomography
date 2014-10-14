//
//  main.cpp
//  tomography
//
//  Created by Даниил Спирин on 29.05.13.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

#include "reconstruction.h"
#include "integration.h"
#include "distribution.h"
#include "dimension.h"
#include "tools.h"
#include "solution.h"

bool useSecondGrid = true;
bool calcParametrs = true;

int timeStart, timeFinish;

string pathToData = "/Users/imaginary/Documents/Science/2012_103_10/";
string pathToProcessedData = pathToData + "tec_processed_model/";

int main (int argc, const char * argv[]) {

    timeStart = 0, timeFinish = 24;
    
    vector<vector<Ray>> data = getData(pathToData.c_str(), timeStart, timeFinish);
    
    Solution solution;
    solution.setLimits(-10.0, 40.0, 30.0, 70.0, timeStart, timeFinish);

    solution.addData(data);
    for (int i = 0; i < 7; ++i) {
        int j = static_cast<int>(pow(static_cast<double>(2), i)+.005);
        solution.addGrid(2*j, (timeFinish-timeStart)*j/2);
    }
    solution.find();
    solution.print();
    
    return 0;
}

