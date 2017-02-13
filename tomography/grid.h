#pragma once

#include <vector>
#include "vector_sparse.h"
#include "dimension.h"

using namespace std;

class Grid : public vector<double> {
public:
    double  left1, right1,
            left2, right2,
            left3, right3;
    size_t  intervals1, intervals2, intervals3;
    Grid(){};
    void set(Dimension latitude, Dimension longitude, Dimension time);
    VectorSparse basis(double x, double y, double z) const;
    
    double operator ()(double x, double y, double z);
    
private:
    void expand();
    void expander(double &left, double &right, size_t &intervals);
};