#ifndef ISOPARAM_H
#define ISOPARAM_H

#include "solver/mesh/mesh2d.h"
#include "solver/finite-element.h"
#include "solver/finite-element1d.h"

double IMapXXi(struct fe*, Elem2D*, double, double);
double IMapYXi(struct fe*, Elem2D*, double, double);
double IMapXEta(struct fe*, Elem2D*, double, double);
double IMapYEta(struct fe*, Elem2D*, double, double);
double IMap1D(struct fe1d*, Elem1D*, double);
double IMapDt1D(struct fe1d*, Elem1D*, double);

double IMapJ(struct fe*, Elem2D*, double, double);
double IMapCyl(struct fe*, Elem2D*, double, double);
double IMapCyl1D(struct fe1d*, Elem1D*, double);

double IEvalLin2Dx(struct fe*, Elem2D*, int, double, double);
double IEvalLin2Dy(struct fe*, Elem2D*, int, double, double);

#endif

