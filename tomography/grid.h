#pragma once

#include <vector>
#include "vector_sparse.h"
#include "dimension.h"

class Grid : public std::vector<double> {
public:
    Grid(){};
    void set(const Dimension &latitude, const Dimension &longitude, const Dimension &time);
    VectorSparse basis(const double x, const double y, const double z) const;
    double operator ()(const double x, const double y, const double z);
private:
    Dimension latitude, longitude, time;
};