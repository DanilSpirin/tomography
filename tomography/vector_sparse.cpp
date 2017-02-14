#include <iostream>
#include "vector_sparse.h"


void VectorSparse::add(float a, int b) {
    phi.push_back(a);
    number.push_back(b);
}
int VectorSparse::getNumber(int i) const {
    return number.at(i);
}
double VectorSparse::getPhi(int i) const {
    return phi.at(i);
}

VectorSparse operator - (const VectorSparse &a) {
    VectorSparse b;
    for (int i = 0; i < a.getSize(); ++i) {
        b.add(-a.getPhi(i), a.getNumber(i));
    }
    return b;
}
VectorSparse operator - (const VectorSparse &a, const VectorSparse &b) {
    if (b.getSize() == 0 || b.getSize() == -1) {
        return a;
    }
    if (a.getSize() == 0 || a.getSize() == -1) {
        return -b;
    }
    VectorSparse c;
    int i1 = 0, i2 = 0;
    while(i1 < a.getSize() || i2 < b.getSize()) {
		if(i1 < a.getSize() && (i2 == b.getSize() || a.number[i1] < b.number[i2])) {
			c.add(a.phi[i1], a.number[i1]);
			i1++;
		}
        else if (i2 < b.getSize() && (i1 == a.getSize() || b.number[i2] < a.number[i1])) {
			c.add(-b.phi[i2], b.number[i2]);
			i2++;
		}
        else {
			c.add(a.phi[i1]-b.phi[i2], a.number[i1]);
			i1++;
			i2++;
		}
    }
    return c;
}

VectorSparse operator * (const double a, const VectorSparse b) {
    VectorSparse tmp(b);
    for (int i = 0; i < tmp.getSize(); ++i) {
        tmp.phi[i] = b.phi[i]*a;
    }
    return tmp;
}

VectorSparse operator * (const VectorSparse b, const double a) {
    VectorSparse tmp(b);
    for (int i = 0; i < tmp.getSize(); ++i) {
        tmp.phi[i] = b.phi[i]*a;
    }
    return tmp;
}