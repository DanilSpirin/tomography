//
//  transformation.cpp
//  model
//
//  Created by Даниил Спирин on 12.09.12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include <cmath>

#include "transformation.h"


void DecartToGeographic::forward(point &R) {
    point r;
    r.R[2] = sqrt(R.R[0] * R.R[0] + R.R[1] * R.R[1] + R.R[2] * R.R[2]);
    r.R[0] = atan2(R.R[1], R.R[0]);
    r.R[1] = asin(R.R[2] / r.R[2]);
    R = r;
}

void DecartToGeographic::backward(point &R) {
    point r;
    r.R[0] = R.R[2] * cos(R.R[1]) * cos(R.R[0]);
    r.R[1] = R.R[2] * cos(R.R[1]) * sin(R.R[0]);
    r.R[2] = R.R[2] * sin(R.R[1]);
    R = r;
}