#pragma once

#include <vector>

#include "dimension.h"
#include "vector_sparse.h"

class Grid : public std::vector<double> {
public:
    void set(const Dimension &latitude, const Dimension &longitude, const Dimension &time);
    VectorSparse basis(const double x, const double y, const double z) const;
    double operator ()(const double x, const double y, const double z) const;
private:
    Dimension latitude, longitude, time;
};
