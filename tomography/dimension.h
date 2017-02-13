#pragma once

struct Dimension {
    double left, right;
    int intervals;
    Dimension();
    Dimension(double left, double right, int itervals, bool inDegrees = true);
    Dimension(Dimension const &foo);
    void toRadian();
    void toDegrees();
private:
    bool inDegrees;
};