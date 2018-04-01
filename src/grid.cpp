#include "grid.hpp"

#include "spline.hpp"

const std::vector<float> spline_y = {0, 0.25, 1, 0.25, 0};
const std::vector<float> spline_x = {-2, -1, 0, 1, 2};
const Spline spline(spline_x, spline_y);

Grid::Grid(const Dimension &latitude, const Dimension &longitude, const Dimension &time)
    : latitude(latitude), longitude(longitude), time(time) {
    this->latitude.expand();
    this->longitude.expand();
    this->latitude.to_radian();
    this->longitude.to_radian();
    data.resize((this->latitude.size() + 1) * (this->longitude.size() + 1) * (this->time.size() + 1), 0);
}

VectorSparse Grid::basis(const float x, const float y, const float z) const {
    VectorSparse basis_vector;

    for (const auto& lat : this->latitude.sequence(x)) {
        for (const auto& lon : this->longitude.sequence(y)) {
            for (const auto& time : this->time.sequence(z)) {
                const float value = spline(this->latitude.get(x, lat))
                                    * spline(this->longitude.get(y, lon))
                                    * spline(this->time.get(z, time));
                const unsigned number = (lat * (this->longitude.size() + 1) + lon) * (this->time.size() + 1) + time;
                basis_vector.push_back(number, value);
            }
        }
    }
    return basis_vector;
}

float Grid::operator()(const float x, const float y, const float z) const {
    float sum = 0;
    VectorSparse temp = basis(x, y, z);
    for (const auto& it : temp) {
        sum += it.value * data[it.index];
    }
    return sum;
}
