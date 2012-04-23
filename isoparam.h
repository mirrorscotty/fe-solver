#ifndef ISOPARAM_H
#define ISOPARAM_H

#include "mesh.h"
#include "finite-element.h"

double IMapXXi(struct fe*, Elem2D*, double, double);
double IMapYXi(struct fe*, Elem2D*, double, double);
double IMapXEta(struct fe*, Elem2D*, double, double);
double IMapYEta(struct fe*, Elem2D*, double, double);

double IMapJ(struct fe*, Elem2D*, double, double);

#endif

