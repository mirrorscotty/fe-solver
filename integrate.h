#ifndef INTEGRATE_H
#define INTEGRATE_H

double quad3(double (*)(double));
double quad300(double (*)(double), double (*)(double));
double quad311(double (*)(double), double (*)(double));
double quad301(double (*)(double), double (*)(double));
double diff(double (*)(double), double);

#endif
