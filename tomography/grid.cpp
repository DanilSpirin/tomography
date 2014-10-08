//
//  grid.cpp
//  vector_3D_task
//
//  Created by Даниил Спирин on 11.05.12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "grid.h"
#include "spline.h"
#include <cmath>
#include "tools.h"

void Grid::set(Dimension latitude, Dimension longitude, Dimension time) {
    clear();
    left1 = latitude.left; right1 = latitude.right; intervals1 = latitude.intervals;
    left2 = longitude.left; right2 = longitude.right; intervals2 = longitude.intervals;
    left3 = time.left; right3 = time.right; intervals3 = time.intervals;
    expand();
    resize((intervals1 + 1) * (intervals2 + 1) * (intervals3 + 1), 0);
}

VectorSparse Grid::basis(double x, double y, double z) const {
    VectorSparse basis_vector;
    
    Spline spline;
    double spline_y[] = {0, 0.25, 1, 0.25, 0};
    double spline_x[] = {-2, -1, 0, 1, 2};
    spline.build(spline_x, spline_y, 5);
    
    int point_i = (int)floor((x-left1)/(right1-left1)*intervals1);
    int point_j = (int)floor((y-left2)/(right2-left2)*intervals2);
    int point_k = (int)floor((z-left3)/(right3-left3)*intervals3);

    int point_i1 = (point_i-1) < 0 ? 0 : point_i - 1;
    int point_i2 = (point_i+2) > intervals1 ? (int)intervals1 : point_i + 2;
    
    int point_j1 = (point_j-1) < 0 ? 0 : point_j - 1;
    int point_j2 = (point_j+2) > intervals2 ? (int)intervals2 : point_j + 2;
    
    int point_k1 = (point_k-1) < 0 ? 0 : point_k - 1;
    int point_k2 = (point_k+2) > intervals3 ? (int)intervals3 : point_k + 2;
    
    for (int i = point_i1; i <= point_i2; ++i) {
        for (int j = point_j1; j <= point_j2; ++j) {
            for (int k = point_k1; k <= point_k2; ++k) {
                float value = 
                  spline((x-(left1+i*(right1-left1)/intervals1))*intervals1/(right1-left1))
                * spline((y-(left2+j*(right2-left2)/intervals2))*intervals2/(right2-left2))
                * spline((z-(left3+k*(right3-left3)/intervals3))*intervals3/(right3-left3));
                unsigned long number = (i * (intervals2 + 1) + j) * (intervals3 + 1) + k;
                basis_vector.add(value, (int)number);
            }
        }
    }
    return basis_vector;
}

double Grid::operator ()(double x, double y, double z) {
    double sum = 0;
    VectorSparse temp = basis(x, y, z);
    for (int j = 0; j < temp.getSize(); ++j) {
        sum += temp.getPhi(j) * operator[](temp.getNumber(j));
    }
    return sum;
}

void Grid::expander(double &left, double &right, size_t &intervals) {
    float step = (right - left) / intervals;
    int extra = 16;
    int m = ceil(extra/step);
    float D = m * step;
    left -= D; right += D;
    
    left = degreeToRadian(left);
    right = degreeToRadian(right);
    intervals += 2*m;
}

void Grid::expand(){
    expander(left1, right1, intervals1);
    expander(left2, right2, intervals2);
}