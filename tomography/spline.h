//
//  spline.h
//  vector_3D_task
//
//  Created by Даниил Спирин on 11.05.12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef vector_3D_task_spline_h
#define vector_3D_task_spline_h

class Spline{
public:
    Spline();
    ~Spline();
    void build(const double *x, const double *y, std::size_t n); // построение сплайна по n точкам (x, y)
    double operator() (double x); // возвращает значение в точке x
private:
    //структура, описывающая сплайн на каждом сегменте сетки
    struct SplineTuple{
        double a, b, c, d, x;
    };
    SplineTuple *splines;
    std::size_t n;
    void freeMemory();
};

#endif