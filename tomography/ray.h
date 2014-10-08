//
//  ray.h
//  vector_3D_task
//
//  Created by Даниил Спирин on 11.05.12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef vector_3D_task_ray_h
#define vector_3D_task_ray_h

#include "point.h"

class BaseRay {
public:
    float integral, time;
    point station, satellite;
};

class Ray : public BaseRay{
public:
    float thetta, phi, angle;
    Ray(const point station, const point satellite, const float time);
    Ray();
    point cross;
    void computeParameters();
private:
    void computeCross();
    void computeAngles();
};

std::istream& operator >> (std::istream&, Ray&);
std::ostream& operator << (std::ostream&, Ray);

#endif