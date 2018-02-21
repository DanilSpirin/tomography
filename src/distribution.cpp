#include <cmath>

#include "distribution.h"
#include "point.h"

Spot::Spot(const point &location, const float peak, const float period, const float intensity, const float size)
    : location(location), peak(peak), period(period), intensity(intensity), size(size) {
    this->location.R[0] /= (180.0f / pi);
    this->location.R[1] /= (180.0f / pi);
}

float Spot::operator()(const point& R, const float time) const {
    return size * std::sin(pi / 2 + (peak - time / 3600) * period) * std::exp(-(R - location).radius_squared() * intensity);
}

Wave::Wave(const point &location, const float start, const float period, const float speed)
    : location(location), start(start), period(period), speed(speed) {
    this->location.R[0] /= (180.0f / pi);
    this->location.R[1] /= (180.0f / pi);
}

float Wave::p(const float r, const float v, const float T) const {
    float r0 = v * T / 2;
    if (r < r0) {
        return 1.0f / 2.0f * (1 - std::cos(pi * r / r0));
    } else {
        return 1;
    }
}

float Wave::f(const float t, const float T, const float n) const {
    return std::sin(2 * pi * t / T) * std::exp(-t * t / 2 / n / n / T / T);
}

float Wave::operator()(const point& R_a, const float time) const {
    point R = R_a;
    point Rc = location;
    DecartToGeographic transformation;
    transformation.backward(R);
    transformation.backward(Rc);
    const float r = (R - Rc).length();
    const float n = 2;
    return p(r, speed, period) * f(r / speed - (time - start), period, n);
}

float ElectronDensityDistribution::operator() (point R, const float t) const {
    if (coordinateTransformation) {
        coordinateTransformation->forward(R);
    }
    return value(R, t);
}

ChepmanLayer::ChepmanLayer() : nmin(0.4f), nm(1.1f), hm(300), H(100), d(102), dt(0) {}

ChepmanLayer::~ChepmanLayer() {}

void ChepmanLayer::add_spot(const point &location, const float peak, const float period, const float intensity,
                            const float size) {
    spots.emplace_back(Spot(location, peak, period, intensity, size));
}

void ChepmanLayer::add_wave(const point &location, const float start, const float period, const float speed) {
    waves.emplace_back(Wave(location, start, period, speed));
}

float ChepmanLayer::value(const point &R, const float time) const {
    const float longitute = R.R[0], latitude = R.R[1], h = R.R[2] - Re;

    const float UT = (time - dt) / 60.0f / 60.0f;    // Universal Time

    const float declination = std::asin(std::sin(23.45f / 180.0f * pi) * std::sin(2 * pi / 365.0f * (d - 82.0f))); // convert to radians
    const float angle = (15 * (UT - 12.0f)) / 180.0f * pi + longitute; // angle - time zone
    const float zenith = std::sin(latitude) * std::sin(declination) + std::cos(latitude) * std::cos(declination) * std::cos(angle);
    const float ksi = (h - hm) / H ;

    const float Q = std::max(nm * zenith, nmin);

    float wavesAddition = 0, spotsAddition = 0;
    for (const auto& i : waves) {
        wavesAddition += i(R, time);
    }

    for (const auto& i : spots) {
        spotsAddition += i(R, time);
    }

    return (Q + spotsAddition + 0.5f * nmin * wavesAddition) * std::exp(1.0f - ksi - std::exp(-ksi)) / 10.0f;
}
