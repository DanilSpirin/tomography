#pragma once

#include "point.h"


class CoordinateTransformation {
public:
    virtual void forward(point &R) = 0;
    virtual void backward(point &R) = 0;
    virtual ~CoordinateTransformation() = default;
};

class DecartToGeographic : public CoordinateTransformation {
public:
    void forward(point &R) override;
    void backward(point &R) override;
    ~DecartToGeographic() override = default;
};
