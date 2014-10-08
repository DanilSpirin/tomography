//
//  integration.h
//  model
//
//  Created by Даниил Спирин on 12.09.12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef model_integration_h
#define model_integration_h

#include "ray.h"
#include "distribution.h"


class Integration {
public:
    virtual double operator()(const BaseRay ray, const ElectronDensityDistribution *model) = 0;
};

class Trapezium : public Integration {
public:
    double operator() (const BaseRay ray, const ElectronDensityDistribution *model);
};

class Rectangle : public Integration {
public:
    double operator() (const BaseRay ray, const ElectronDensityDistribution *model);
};
#endif
