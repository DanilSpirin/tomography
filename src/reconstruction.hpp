#pragma once

#include <vector>

#include "grid.hpp"

class Solver {
public:
    Solver() = default;

    virtual void operator()(Grid &grid, const std::vector<VectorSparse> &a, const std::vector<float> &m, bool onlyPositive=true) const = 0;
    virtual ~Solver() = default;
};


class Art : public Solver {
public:
    void operator()(Grid &grid, const std::vector<VectorSparse> &a, const std::vector<float> &m, bool onlyPositive=true) const;
};

class Sirt : public Solver {
public:
    void operator()(Grid &grid, const std::vector<VectorSparse> &a, const std::vector<float> &m, bool onlyPositive=true) const;
};
