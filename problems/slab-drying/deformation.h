#ifndef DEFORMATION_H
#define DEFORMATION_H

#include "finite-element1d.h"

double DeformationGrad(struct fe1d*, double, double);
double DeformGradPc(struct fe1d*, double, double);

#endif

