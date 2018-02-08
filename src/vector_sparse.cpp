#include <iostream>

#include "vector_sparse.h"

size_t VectorSparse::size() const {
    if (_value.size() == _index.size()) {
        return _value.size();
    } else {
        throw "Problem in vector sparse. Index and value sizes differ.";
    }
}

void VectorSparse::push_back(const unsigned index, const double value) {
    _index.push_back(index);
    _value.push_back(value);
}

void VectorSparse::push_back(const Element &element) {
    _index.push_back(element.index);
    _value.push_back(element.value);
}

Element VectorSparse::operator[](const unsigned i) const {
    return Element(_index[i], _value[i]);
}

Element VectorSparse::at(const unsigned i) const {
    return Element(_index.at(i), _value.at(i));
}

VectorSparse operator - (const VectorSparse &a) {
    VectorSparse result(a);
    for (unsigned i = 0; i < result.size(); ++i) {
        result._value[i] = -result._value[i];
    }
    return result;
}

VectorSparse operator - (const VectorSparse &a, const VectorSparse &b) {
    const auto left_size = a.size();
    const auto right_size = b.size();
    if (!right_size) {
        return a;
    }
    if (!left_size) {
        return -b;
    }
    VectorSparse result;
    unsigned i_left = 0;
    unsigned i_right = 0;
    while (i_left < left_size && i_right < right_size) {
        const auto left = a[i_left];
        const auto right = b[i_right];
        if (left.index < right.index) {
            result.push_back(left);
            ++i_left;
        } else if (right.index < left.index) {
            result.push_back(-right);
            ++i_right;
        } else {
            result.push_back(left.index, left.value - right.value);
            ++i_left;
            ++i_right;
        }
    }
    while (i_left < left_size) {
        const auto left = a[i_left];
        result.push_back(left);
        ++i_left;
    } 
    while (i_right < right_size) {
        const auto right = b[i_right];
        result.push_back(-right);
        ++i_right;
    }
    return result;
}

VectorSparse operator * (const double a, const VectorSparse &b) {
    VectorSparse result(b);
    for (unsigned i = 0; i < result.size(); ++i) {
        result._value[i] *= a;
    }
    return result;
}

VectorSparse operator * (const VectorSparse &a, const double b) {
    VectorSparse result(a);
    for (unsigned i = 0; i < result.size(); ++i) {
        result._value[i] *= b;
    }
    return result;
}

VectorSparse operator / (const VectorSparse &a, const double b) {
    VectorSparse result(a);
    for (unsigned i = 0; i < result.size(); ++i) {
        result._value[i] /= b;
    }
    return result;
}
