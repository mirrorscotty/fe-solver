#ifndef HEAT_TRANSFER_H
#define HEAT_TRANSFER_H

#include "matrix.h"
#include "mesh1d.h"

struct fe1d;

double ResHeat(struct fe1d *, matrix *, Elem1D *, double, int, int);
double ResDtHeat(struct fe1d *, matrix *, Elem1D *, double, int, int);
double ConvBCHeat(struct fe1d *, int);
double ExternalTemp(struct fe1d *, int);

#endif

