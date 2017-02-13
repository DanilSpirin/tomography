#include "dimension.h"
#include "ray.h"

Dimension::Dimension(double left, double right, int intervals, bool isInDegrees) {
    this->left = left;
    this->right = right;
    this->intervals = intervals;
    this->inDegrees = isInDegrees;
}

Dimension::Dimension(Dimension const &foo) {
    this->left = foo.left;
    this->right = foo.right;
    this->intervals = foo.intervals;
    this->inDegrees = foo.inDegrees;
}

Dimension::Dimension() {
    
}

void Dimension::toRadian() {
    if (inDegrees) {
        left = left / 180 * pi;
        right = right / 180 * pi;
        inDegrees = false;
    }
}

void Dimension::toDegrees() {
    if (!inDegrees) {
        left = left / pi * 180;
        right = right / pi * 180;
        inDegrees = true;
    }
}