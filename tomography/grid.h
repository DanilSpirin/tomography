#pragma once

#include <vector>
#include "vector_sparse.h"
#include "dimension.h"

using namespace std;

class Grid : public vector<double> {
public:
    Grid(){};
    void set(const Dimension &latitude, const Dimension &longitude, const Dimension &time);
    VectorSparse basis(const double x, const double y, const double z) const;
    double operator ()(const double x, const double y, const double z);
private:
    Dimension latitude, longitude, time;
};