#pragma once

#include <vector>

struct Element {
    unsigned index;
    double value;
    Element(const unsigned index, const double value) : index(index), value(value) {}
    Element operator-() const {return Element(index, -value);}
};

class VectorSparse {
    std::vector<double> _value;
    std::vector<unsigned> _index;
public:
    size_t size() const;
    void push_back(const unsigned index, const double value);
    void push_back(const Element &element);
    Element at(const unsigned i) const;
    Element operator[](const unsigned i) const;
    friend VectorSparse operator - (const VectorSparse& a, const VectorSparse& b);
    friend VectorSparse operator - (const VectorSparse& a);
    friend VectorSparse operator * (const double a, const VectorSparse &b);
    friend VectorSparse operator * (const VectorSparse &b, const double a);
    friend VectorSparse operator / (const VectorSparse &a, const double b);
};