#include <cmath>

#include "distribution.h"
#include "ray.h"

Spot::Spot(const point &location, const double peak, const double period, const double intensity, const double size)
    : location(location), peak(peak), period(period), intensity(intensity), size(size) {
    this->location.R[0] /= (180.0 / pi);
    this->location.R[1] /= (180.0 / pi);
}

double Spot::operator()(const point& R, const double time) const {
    return size * sin(pi / 2 + (peak - time / 3600) * period) * exp(-(R - location).radius_squared() * intensity);
}

Wave::Wave(const point &location, const double start, const double period, const double speed)
    : location(location), start(start), period(period), speed(speed) {
    this->location.R[0] /= (180.0 / pi);
    this->location.R[1] /= (180.0 / pi);
}

double Wave::p(const double r, const double v, const double T) const {
    double r0 = v * T / 2;
    if (r < r0) {
        return 1 / 2 * (1 - cos(pi * r / r0));
    } else {
        return 1;
    }
}

double Wave::f(const double t, const double T, const double n) const {
    return sin(2 * pi * t / T) * exp(-t * t / 2 / n / n / T / T);
}

double Wave::operator()(const point& R_a, const double time) const {
    point R = R_a;
    point Rc = location;
    DecartToGeographic transformation;
    transformation.backward(R);
    transformation.backward(Rc);
    const double r = (R - Rc).length();
    const double n = 2;
    return p(r, speed, period) * f(r / speed - (time - start), period, n);
}
double ElectronDensityDistribution::operator() (point R, const double t) const {
    if (coordinateTransformation) {
        coordinateTransformation->forward(R);
    }
    return value(R, t);
}

ChepmanLayer::ChepmanLayer() : nmin(0.4), nm(1.1), hm(300), H(100), d(102), dt(0) {}

ChepmanLayer::~ChepmanLayer() {}

void ChepmanLayer::add_spot(const point &location, const double peak, const double period, const double intensity,
                            const double size) {
    spots.push_back(Spot(location, peak, period, intensity, size));
}

void ChepmanLayer::add_wave(const point &location, const double start, const double period, const double speed) {
    waves.push_back(Wave(location, start, period, speed));
}

double ChepmanLayer::value(const point &R, const double time) const {
    const double longitute = R.R[0], latitude = R.R[1], h = R.R[2] - Re;

    const double UT = (time - dt) / 60.0 / 60.0;    //всемирное время

    const double declination = asin(sin(23.45 / 180.0 * pi) * sin(2 * pi / 365.0 * (d - 82.0))); //осуществлен перевод в радианы
    const double angle = (15 * (UT - 12.0)) / 180.0 * pi + longitute; //осуществлен перевод в радианы. angle - ЧАСОВОЙ ПОЯС
    const double zenith = sin(latitude) * sin(declination) + cos(latitude) * cos(declination) * cos(angle);
    const double ksi = (h - hm) / H ;

    const double Q = std::max(nm * zenith, nmin);

    double wavesAddition = 0, spotsAddition = 0;
    for (const auto& i : waves) {
        wavesAddition += i(R, time);
    }

    for (const auto& i : spots) {
        spotsAddition += i(R, time);
    }

    return (Q + spotsAddition + 0.5 * nmin * wavesAddition) * exp(1 - ksi - exp(-ksi)) / 10;
}
