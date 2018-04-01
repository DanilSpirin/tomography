#pragma once

#include "dimension.hpp"
#include "vector_sparse.hpp"

class Grid {
public:
    Grid(const Dimension &latitude, const Dimension &longitude, const Dimension &time);
    VectorSparse basis(const float x, const float y, const float z) const;

    std::size_t size() const { return data.size(); }

    float operator()(const float x, const float y, const float z) const;

    float& operator[](const std::size_t i) { return  data[i]; }
    const float& operator[](const std::size_t i) const { return data[i]; }

    auto begin() {return data.begin();}
    auto end() { return data.end(); }
    auto begin() const { return data.begin(); }
    auto end() const { return data.end(); }
    auto cbegin() const { return data.cbegin(); }
    auto cend() const { return data.cend(); }
private:
    std::vector<float> data;
    Dimension latitude, longitude, time;
};
