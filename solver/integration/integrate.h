#ifndef INTEGRATE_H
#define INTEGRATE_H

#include "solver/integration/basis.h"
#include "solver/mesh/mesh2d.h"
#include "solver/finite-element.h"
#include "solver/finite-element1d.h"


double quad2d3generic(struct fe*, matrix*, Elem2D*,
                      double (*)(struct fe*, matrix*, Elem2D*, double, double, int, int),
                      int, int);
double quad1d3generic(struct fe1d*, matrix*, Elem1D*,
                      double (*)(struct fe1d*, matrix*, Elem1D*, double, int, int),
                      int, int);

/* Depreciated */
double quad3(double (*)(double));
double quad300(double (*)(double), double (*)(double));
double quad2d3(struct fe*, Elem2D*, int, int, int);
double quad2d32d3(struct fe*, Elem2D*, int, int, int, int);

#endif

