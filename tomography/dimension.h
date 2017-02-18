#pragma once

struct Dimension {
    double left, right;
    Dimension();
    Dimension(const double left, const double right, const size_t itervals, const bool inDegrees = true);
    Dimension(Dimension const &foo);
    void toRadian();
    void toDegrees();
    void expand(const uint16_t extra = 16);
    size_t size() const;
    std::vector<size_t> sequence(const double i) const; // TODO: come up withe a better name
    double get(const double x, const double i) const; // TODO: come up with a better name
    friend std::ostream& operator << (std::ostream& out, const Dimension& dimension);
private:
    bool inDegrees;
    double step;
    size_t intervals;
};