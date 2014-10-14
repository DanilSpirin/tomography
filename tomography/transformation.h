//
//  transformation.h
//  model
//
//  Created by Даниил Спирин on 12.09.12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef model_transformation_h
#define model_transformation_h

#include "point.h"


class CoordinateTransformation {
public:
    virtual void forward(point &R) = 0;
    virtual void backward(point &R) = 0;
    virtual ~CoordinateTransformation();
};

class DecartToGeographic : public CoordinateTransformation {
public:
    void forward(point &R);
    void backward(point &R);
};

#endif
