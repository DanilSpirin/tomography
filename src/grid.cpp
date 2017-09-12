#include <cmath>
#include <numeric>

#include "grid.h"
#include "spline.h"
#include "tools.h"

void Grid::set(const Dimension &latitude, const Dimension &longitude, const Dimension &time) {
    clear();
    this->latitude = latitude;
    this->longitude = longitude;
    this->time = time;
    this->latitude.expand();
    this->longitude.expand();
    this->latitude.toRadian();
    this->longitude.toRadian();
    resize((this->latitude.size() + 1) * (this->longitude.size() + 1) * (this->time.size() + 1), 0);
}

VectorSparse Grid::basis(const double x, const double y, const double z) const {
    VectorSparse basis_vector;

    Spline spline;
    const double spline_y[] = {0, 0.25, 1, 0.25, 0};
    const double spline_x[] = {-2, -1, 0, 1, 2};
    spline.build(spline_x, spline_y, 5);

    for (const auto& lat : this->latitude.sequence(x)) {
        for (const auto& lon : this->longitude.sequence(y)) {
            for (const auto& time : this->time.sequence(z)) {
                const double value = spline(this->latitude.get(x, lat))
                                    * spline(this->longitude.get(y, lon))
                                    * spline(this->time.get(z, time));
                const unsigned number = (lat * (this->longitude.size() + 1) + lon) * (this->time.size() + 1) + time;
                basis_vector.push_back(number, value);
            }
        }
    }
    return basis_vector;
}

double Grid::operator()(const double x, const double y, const double z) {
    double sum = 0;
    VectorSparse temp = basis(x, y, z);
    for (unsigned i = 0; i < temp.size(); ++i) {
        const auto element = temp[i];
        sum += element.value * operator[](element.index);
    }
    return sum;
}
