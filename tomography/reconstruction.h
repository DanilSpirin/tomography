//
//  rt.h
//  vector_3D_task
//
//  Created by Даниил Спирин on 11.05.12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef vector_3D_task_rt_h
#define vector_3D_task_rt_h

#include <vector>
#include "vector_sparse.h"
#include "grid.h"
using namespace std;

void iterationArt(Grid &grid, const vector<VectorSparse> &a, const vector<float> &m, bool onlyPositive=true);
void iterationSirt(Grid &grid, const vector<VectorSparse> &a, const vector<float> &m, bool onlyPositive=true);

#endif