#pragma once

#include <vector>
using namespace std;

class VectorSparse {
    vector<float> phi;
    vector<int> number;
public:
    size_t getSize() const {return phi.size() == number.size() ? phi.size() : -1;};
    void add(float a, int b);
    double getPhi(int i) const;
    int getNumber(int i) const;
    VectorSparse& operator *= (const double a);
    friend VectorSparse operator - (const VectorSparse& a, const VectorSparse& b);
    friend VectorSparse operator - (const VectorSparse& a);
    friend VectorSparse operator * (const double a, const VectorSparse b);
    friend VectorSparse operator * (const VectorSparse b, const double a);
};