#ifndef INTEGRATE_H
#define INTEGRATE_H

#include "basis.h"
#include "mesh.h"
#include "finite-element.h"

double quad3(double (*)(double));
double quad300(double (*)(double), double (*)(double));
double quad311(double (*)(double), double (*)(double));
double quad301(double (*)(double), double (*)(double));
double diff(double (*)(double), double);

double diff2dx(basis*, int, double, double);
double diff2dy(basis*, int, double, double);
double quad2d3(struct fe*, Elem2D*, int, int, int);
double quad2d32d3(struct fe*, Elem2D*, int, int, int, int);
double quad2d3generic(struct fe*, matrix*, Elem2D*,
                      double (*)(struct fe*, matrix*, Elem2D*, double, double, int, int),
                      int, int);

#endif
