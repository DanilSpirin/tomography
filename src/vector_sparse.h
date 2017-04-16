#pragma once

#include <vector>

class VectorSparse {
    std::vector<double> phi;
    std::vector<int> number;
public:
    size_t getSize() const {return phi.size() == number.size() ? phi.size() : -1;};
    void add(double a, int b);
    double getPhi(int i) const;
    int getNumber(int i) const;
//    VectorSparse& operator *= (const double a);
    friend VectorSparse operator - (const VectorSparse& a, const VectorSparse& b);
    friend VectorSparse operator - (const VectorSparse& a);
    friend VectorSparse operator * (const double a, const VectorSparse b);
    friend VectorSparse operator * (const VectorSparse b, const double a);
    friend VectorSparse operator / (const VectorSparse a, const double b);
};