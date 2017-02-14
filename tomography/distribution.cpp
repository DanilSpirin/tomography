#include <iostream>
#include <cmath>

#include "distribution.h"
#include "ray.h"

using namespace std;

Spot::Spot() {}
Spot::~Spot() {}

Spot::Spot(point location, double peak, double period, double intensity, double size) {
    this->location = location;
    this->location.R[0] /= (180.0/pi);
    this->location.R[1] /= (180.0/pi);
    this->peak = peak;
    this->period = period;
    this->intensity = intensity;
    this->size = size;
}

Spot::Spot(const Spot &a)  {
    this->location = a.location;
    this->peak = a.peak;
    this->period = a.period;
    this->intensity = a.intensity;
    this->size = a.size;
}

double Spot::operator()(const point R, const double time) const {
    return size*sin(pi/2+(peak-time/3600)*period)*exp((-(R-location).R[0]*(R-location).R[0] - (R-location).R[1]*(R-location).R[1])*intensity);
}

Wave::Wave() {}
Wave::~Wave() {}

Wave::Wave(point location, double start, double period, double speed) {
    this->location = location;
    this->location.R[0] /= (180.0/pi);
    this->location.R[1] /= (180.0/pi);
    this->start = start;
    this->period = period;
    this->speed = speed;
}

Wave::Wave(const Wave &a) {
    this->location = a.location;
    this->start = a.start;
    this->period = a.period;
    this->speed = a.speed;
}
double Wave::p(const double r, const double v, const double T) const {
    double r0 = v * T / 2;
    if (r < r0) {
        return 1 / 2 * (1 - cos(pi * r / r0));
    }
    else {
        return 1;
    }
}

double Wave::f(const double t, const double T, const double n) const {
    return sin(2*pi*t/T)*exp(-t*t/2/n/n/T/T);
}

double Wave::operator()(const point R_a, const double time) const {
    point R = R_a;
    point Rc = location;
    DecartToGeographic transformation;
    transformation.backward(R);
    transformation.backward(Rc);
    double r = sqrt((R.R[0]-Rc.R[0])*(R.R[0]-Rc.R[0]) + (R.R[1]-Rc.R[1])*(R.R[1]-Rc.R[1]) + (R.R[2]-Rc.R[2])*(R.R[2]-Rc.R[2]));
    double n = 2;
    return p(r, speed, period) * f(r/speed - (time-start), period, n); 
}
double ElectronDensityDistribution::operator() (point R, const double t) const {
    if (coordinateTransformation) {
        coordinateTransformation->forward(R);
    }
    return value(R, t);
}

ChepmanLayer::ChepmanLayer() : nmin(0.4), nm(1.1), hm(300), H(100), d(102), dt(0), waves(NULL), spots(NULL) {}

ChepmanLayer::~ChepmanLayer() {}

void ChepmanLayer::addSpot(point location, double peak, double period, double intensity, double size) {
    spots.push_back(Spot(location, peak, period, intensity, size));
}

void ChepmanLayer::addWave(point location, double start, double period, double speed) {
    waves.push_back(Wave(location, start, period, speed));
}

double ChepmanLayer::value(const point R, const double time) const {
    double longitute = R.R[0], latitude = R.R[1], h = R.R[2]-Re;
	
	double zenith;	//zenith = cos зенитного угла
	double angle;	//часовой угол Солнца (на гринвиче) т.е. к нему нужно добавить долготу !
	double declination; //склонение Солнца
	double UT = ((double)time-dt)/60.0/60.0;    //всемирное время
	double ksi;
    
	declination = asin(sin(23.45/180.0*pi)*sin(2*pi/365.0*(double(d) - 82.0))); //осуществлен перевод в радианы
	angle = (15*(UT - 12.0))/180.0*pi + longitute; //осуществлен перевод в радианы. angle - ЧАСОВОЙ ПОЯС 
	zenith = sin(latitude)*sin(declination) + cos(latitude)*cos(declination)*cos(angle);
	ksi = (h - hm)/H ;
    
    double Q = nm * zenith;
	if (Q < nmin) {
        Q = nmin;
    }
    
    double wavesAddition = 0, spotsAddition = 0;
    if (waves.size() > 0) {
        for (auto i = waves.begin(); i != waves.end(); ++i) {
            wavesAddition += (*i)(R, time);
        }
    }

    if (spots.size() > 0) {
        for (auto i = spots.begin(); i != spots.end(); ++i) {
            spotsAddition += (*i)(R, time);
        }
    }
    
    return (Q + spotsAddition + 0.5*nmin*wavesAddition)*exp(1-ksi-exp(-ksi))/10;
}