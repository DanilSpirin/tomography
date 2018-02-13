#pragma once

#include <vector>

#include "grid.h"
#include "vector_sparse.h"

class Solver {
public:
    Solver() {}
    virtual void operator()(Grid &grid, const std::vector<VectorSparse> &a, const std::vector<float> &m, bool onlyPositive=true) const = 0;
    virtual ~Solver() {}
};


class Art : public Solver {
public:
    void operator()(Grid &grid, const std::vector<VectorSparse> &a, const std::vector<float> &m, bool onlyPositive=true) const;
};

class Sirt : public Solver {
public:
    void operator()(Grid &grid, const std::vector<VectorSparse> &a, const std::vector<float> &m, bool onlyPositive=true) const;
};
