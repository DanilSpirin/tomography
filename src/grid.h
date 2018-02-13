#pragma once

#include <vector>

#include "dimension.h"
#include "vector_sparse.h"

class Grid : public std::vector<float> {
public:
    void set(const Dimension &latitude, const Dimension &longitude, const Dimension &time);
    VectorSparse basis(const float x, const float y, const float z) const;
    float operator ()(const float x, const float y, const float z) const;
private:
    Dimension latitude, longitude, time;
};
