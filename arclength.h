#ifndef ARC_LENGTH_H
#define ARC_LENGTH_H

#include "matrix.h"
#include "finite-element.h"

double ALength(struct fe*, matrix*, matrix*, Elem2D*, double, double, double, double);
double ALengthx(struct fe*, matrix*, matrix*, Elem2D*, double, double, double, double);
double ALengthy(struct fe*, matrix*, matrix*, Elem2D*, double, double, double, double);

matrix *AssembleBottomRow(struct fe*, matrix*, matrix*);
matrix *AssembleRightCol(struct fe*, matrix*, matrix*);
double AssembleBotCorner(struct fe*, matrix*, matrix*);

matrix *AddArcLengthConstraint(struct fe*, matrix*, matrix*);

#endif

