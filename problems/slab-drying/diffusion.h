#ifndef DIFFUSION_H
#define DIFFUSION_H

#include "matrix.h"
#include "mesh1d.h"

struct fe1d;

double ResMass(struct fe1d *, matrix *, Elem1D *, double, int, int);
double ResDtMass(struct fe1d *, matrix *, Elem1D *, double, int, int);
double ConvBCMass(struct fe1d *, int);
double ExternalConc(struct fe1d *, int);

#endif

