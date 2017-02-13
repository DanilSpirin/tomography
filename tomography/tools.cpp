#include <iostream>
#include <fstream>
#include <string>

#include <sys/types.h>
#include <dirent.h>

#include "tools.h"
#include "math.h"
#include "ray.h"
#include "distribution.h"
#include "integration.h"

extern int timeStart, timeFinish;
extern string pathToProcessedData;

void dataToSle(const vector<vector<Ray> > &data, vector<VectorSparse> &phi, vector<float> &integrals, const Grid &test) {
    phi.clear();
    integrals.clear();
    for (int j = 0; j < data.size(); ++j) {
        for (int i = 0; i < data[j].size()-1; ++i) {
//          Разностная схема
            phi.push_back((1/cos(data[j][i+1].angle)) * test.basis(data[j][i+1].phi, data[j][i+1].thetta, data[j][i+1].time) - (1/cos(data[j][i].angle)) * test.basis(data[j][i].phi, data[j][i].thetta, data[j][i].time));
            integrals.push_back(data[j][i+1].integral - data[j][i].integral);
        }
    }
}

void dataToSle(const vector<vector<Ray> > &data, vector<VectorSparse> &phi, const Grid &test) {
    phi.clear();
    for (int j = 0; j < data.size(); ++j) {
        for (int i = 0; i < data[j].size()-1; ++i) {
            phi.push_back((1/cos(data[j][i+1].angle)) * test.basis(data[j][i+1].phi, data[j][i+1].thetta, data[j][i+1].time) - (1/cos(data[j][i].angle)) * test.basis(data[j][i].phi, data[j][i].thetta, data[j][i].time));
        }
    }
}

double computeResidual(const Grid &x, const vector<VectorSparse> &A, const vector<float> &m) {
    vector<double> difference(m.size(), 0);
    for (int i = 0; i < A.size(); ++i) {
        for (int j = 0; j < A[i].getSize(); ++j) {
            difference[i] += A[i].getPhi(j)*x[A[i].getNumber(j)];
        }
        difference[i] -= m[i];
    }
    double sum = 0;
    for (int i = 0; i < difference.size(); ++i) {
        sum += difference[i]*difference[i];
    }
    
    return sqrt(sum);
}

vector<float> computeVectorResidual(const Grid &x, const vector<VectorSparse> &A, const vector<float> &m) {
    vector<float> difference(m.size(), 0);
    for (int i = 0; i < A.size(); ++i) {
        for (int j = 0; j < A[i].getSize(); ++j) {
            difference[i] += A[i].getPhi(j) * x[A[i].getNumber(j)];
        }
        difference[i] -= m[i];
    }
    for (int i = 0; i < difference.size(); ++i) {
        difference[i] = -difference[i];
    }
    return difference;
}

// отсеиваем скрытые файлы и файлы не являющиеся данными
int sel(const struct dirent *d) {
    if (!strcmp(".", d->d_name) || !strcmp("..", d->d_name) || !strcmp(".DS_Store", d->d_name)) return 0;
    char* p = strstr(d->d_name, ".dat");
    return p ? 1 : 0;
}
vector<vector<Ray> > getData(const char *pathToData, int startTime, int finishTime) {
    ChepmanLayer chepmanLayer;
    chepmanLayer.coordinateTransformation = new DecartToGeographic;
    startTime *= 3600;
    finishTime *= 3600;
    dirent **nameList;
    int numberOfBundles = scandir(pathToData, &nameList, sel, 0);
    char temp[100];
    vector<vector<Ray> > data;
    
    if (numberOfBundles < 0) {
        perror("Scandir");
    }
    else if (numberOfBundles == 0) {
        cout << "No data files found!" << endl;
    }
    else {
        while (numberOfBundles--) {
            vector<Ray> tempBundle;
            sprintf(temp, "%s%s", pathToData, nameList[numberOfBundles]->d_name);
            ifstream gps(temp);
            
            if (!gps) {
                cout << "Can't open file " << temp << endl;
            }
            else {
                int numberOfRays;
                gps >> numberOfRays;
                
                while (numberOfRays--) {
                    Ray tempRay;
                    gps >> tempRay;
                        if (tempRay.time >= startTime && tempRay.time <= finishTime) {
                            tempRay.computeParameters();
                            tempBundle.push_back(tempRay);
                        }
                }
                gps.close();
            }
            if (tempBundle.size() > 1) {
                data.push_back(tempBundle);
            }
            free(nameList[numberOfBundles]);
        }
    }
    
    free(nameList);
    
    return data;
}

list<pair<double, double>> getStationList(vector<vector<Ray>> data) {
    list<pair<double, double>> stations;
    bool found;
    ChepmanLayer chepmanLayer;
    chepmanLayer.coordinateTransformation = new DecartToGeographic;
    for (auto i = data.begin(); i != data.end(); ++i) {
        for (auto j = (*i).begin(); j != (*i).end(); ++j) {
            point station = (*j).station;
            chepmanLayer.coordinateTransformation->forward(station);
            found = false;
            for (auto it = stations.begin(); it != stations.end(); ++it) {
                if ((*it).first == station.R[0] && (*it).second == station.R[1]) {
                    found = true;
                }
            }
            if (!found) {
                stations.push_back(make_pair(station.R[0], station.R[1]));
            }
        }
    }
    return stations;
}


void solveSle(Grid &grid, const vector<VectorSparse> &matrix, const vector<float> integrals, double error, bool onlyPositive) {
    double initialResidual = computeResidual(grid, matrix, integrals);
    double firstRes, secondRes, currentRes;
    double limit = 0;
    double iterations = 50;
    double counter = 0;
    for (int i = 0; i < iterations; ++i) {
        iterationSirt(grid, matrix, integrals, onlyPositive);
    }
    firstRes = computeResidual(grid, matrix, integrals)/initialResidual;
    for (int i = 0; i < iterations; ++i) {
        iterationSirt(grid, matrix, integrals, onlyPositive);
    }
    secondRes = computeResidual(grid, matrix, integrals)/initialResidual;
    limit = (iterations*2*secondRes - iterations*firstRes)/iterations;
    cout << limit << endl;
    currentRes = secondRes;
    while (currentRes / limit > 1 + error) {
        iterationSirt(grid, matrix, integrals, onlyPositive);
        currentRes = computeResidual(grid, matrix, integrals)/initialResidual;
        counter++;
        if (counter > 500) {
            cout << "stopped at current/limit = " << currentRes/limit << endl;
            break;
        }
        cout << counter << endl;
    }
}


double degreeToRadian(double degree){
    return degree/180*pi;
}

double radianToDegree(double radian) {
    return radian*180/pi;
}

void computeParametrs(Grid &crude, Grid &accurate, vector<VectorSparse> sleMatrix, vector<float> integrals, bool useSecondGrid, ElectronDensityDistribution *model, Dimension latitude, Dimension longitude, Dimension time, int intervals, int intervalsTime, double initialResidual) {
    latitude.toRadian();
    longitude.toRadian();
    double reconstructionSum = 0, modelSum = 0;
    int density = 100; // Количество точек по оси, по которым строится область
    for (int t = timeStart; t < timeFinish; ++t) {
        for (int i = 0; i <= density; ++i) {
            for (int j = 0; j <= density; ++j) {
                double phi = latitude.left + (latitude.right - latitude.left) / density * i;
                double theta = longitude.left + (longitude.right - longitude.left) / density * j;
                double time = t * 3600;
                
                point station(phi, theta, Re), satellite(phi, theta, Re+1000);
                model->coordinateTransformation->backward(station);
                model->coordinateTransformation->backward(satellite);
                
                Ray L(station, satellite, time);
                
                double crudeValue, accurateValue, sumValue, modelValue;
                Rectangle integral;
                
                crudeValue = crude(phi, theta, time);
                modelValue = integral(L, model);
                
                modelSum += modelValue * modelValue;
                if (useSecondGrid) {
                    accurateValue = accurate(phi, theta, time);
                    sumValue = crudeValue + accurateValue;
                }
                if (useSecondGrid) {
                    reconstructionSum += (sumValue - modelValue)*(sumValue - modelValue);
                }
                else {
                    reconstructionSum += (crudeValue - modelValue)*(crudeValue - modelValue);
                }
            }
        }
    }
    
    latitude.toDegrees();
    longitude.toDegrees();

    ofstream parametrs((pathToProcessedData+"parametrs.txt").c_str(), ios::app);
    
    parametrs << intervals << '\t' << intervalsTime << '\t';
    if (useSecondGrid) {
        parametrs << computeResidual(accurate, sleMatrix, integrals)/initialResidual << '\t';
    }
    else {
        parametrs << computeResidual(crude, sleMatrix, integrals)/initialResidual << '\t';
    }
    parametrs << sqrt(reconstructionSum/modelSum) << '\t' << (longitude.right-longitude.left)/intervals << '\t' << (time.right-time.left)/intervalsTime/60<< endl;
    parametrs.close();

}

list<int> createListOfIntervals(int first, int last) {
    list<int> foo;
    do {
        foo.push_back(first);
        first *= 2;
    } while (first <= last);
    return foo;
}
