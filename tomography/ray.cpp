//
//  ray.cpp
//  vector_3D_task
//
//  Created by Даниил Спирин on 11.05.12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <math.h>
#include "ray.h"
#include "transformation.h"

using namespace std;

std::istream& operator >> (std::istream& in, Ray &c) {
    in >> c.station >> c.satellite >> c.time >> c.integral;
    return in;
}
std::ostream& operator << (std::ostream& out, Ray c) {
    out << c.station << endl << c.satellite << " " << c.time << " " << c.integral;
    return out;
}

Ray::Ray() {}

Ray::Ray(point station, point satellite, float time) {
    this->station = station;
    this->satellite = satellite;
    this->time = time;
}
void Ray::computeCross() {
    point dr = satellite - station;
    double rdr = station * dr;
    double t = (-rdr + sqrt(rdr * rdr + (dr * dr) * ((Re + h) * (Re + h) - (station * station)))) / (dr * dr);
	cross = station + t * (satellite - station);
}
void Ray::computeAngles() {
    double ac = (cross * (satellite - station)) / sqrt((cross * cross) * ((satellite - station) * (satellite - station)));
    if (ac > 1) ac = 1;
    angle = acos(ac);
    thetta = asin(cross.R[2] / (Re + h));
    phi = atan2(cross.R[1], cross.R[0]);
}
void Ray::computeParameters() {
    computeCross();
    computeAngles();
}
