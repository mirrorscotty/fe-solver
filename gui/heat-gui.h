#ifndef HEAT_GUI_H
#define HEAT_GUI_H

#include "matrix.h"
#include "mesh1d.h"

struct fe1d;

double Residual(struct fe1d *, matrix *, Elem1D *, double, int, int);
double ResDt(struct fe1d *, matrix *, Elem1D *, double, int, int);
matrix* CreateElementMatrix(struct fe1d *, Elem1D *, matrix *);
matrix* CreateDTimeMatrix(struct fe1d *, Elem1D *, matrix *);
matrix* CreateElementLoad(struct fe1d *, Elem1D *, matrix *);
int IsOnRightBoundary(struct fe1d *, int);
int IsOnLeftBoundary(struct fe1d *, int);
double Left(struct fe1d *, int);
double Right(struct fe1d *, int);
double Zero(struct fe1d *, int);
double ConvBC(struct fe1d *, int);
void ApplyAllBCs(struct fe1d *);
double InitTemp(double);
double InitC(double);
double react1(double, double, double);
double react2(double, double, double);


void seth(double);
double printh(void);

#endif

