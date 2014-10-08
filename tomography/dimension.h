//
//  dimension.h
//  model
//
//  Created by Даниил Спирин on 29.05.13.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#ifndef model_dimension_h
#define model_dimension_h

struct Dimension {
    double left, right;
    int intervals;
    Dimension();
    Dimension(double left, double right, int itervals, bool inDegrees = true);
    Dimension(Dimension const &foo);
    void toRadian();
    void toDegrees();
private:
    bool inDegrees;
};

#endif
