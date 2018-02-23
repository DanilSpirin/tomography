#include <iostream>

#include "vector_sparse.hpp"


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
