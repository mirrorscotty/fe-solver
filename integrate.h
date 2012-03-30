#ifndef INTEGRATE_H
#define INTEGRATE_H

#include "basis.h"

double quad3(double (*)(double));
double quad300(double (*)(double), double (*)(double));
double quad311(double (*)(double), double (*)(double));
double quad301(double (*)(double), double (*)(double));
double diff(double (*)(double), double);

double diff2dx(basis*, int, double, double);
double diff2dy(basis*, int, double, double);
double quad2d3(basis*, int, int, int);

#endif
