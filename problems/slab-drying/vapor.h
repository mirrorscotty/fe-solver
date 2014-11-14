#ifndef DIFFUSION_H
#define DIFFUSION_H

#include "matrix.h"
#include "mesh1d.h"

struct fe1d;

double ResVap(struct fe1d *, matrix *, Elem1D *, double, int, int);
double ResDtVap(struct fe1d *, matrix *, Elem1D *, double, int, int);
double ConvBCVap(struct fe1d *, int);
double ExternalConcVap(struct fe1d *, int);

#endif

