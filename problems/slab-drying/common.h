#ifndef COMMON_H
#define COMMON_H

#define TVAR 0
#define CVAR 1 
#define PVAR 1

#define TAMB 320 // K
#define TINIT 298 //K
#define HCONV 500

#define CAMB .05 // kg/kg db
#define CINIT .2 // kg/kg db
#define KC_CONV 2e-8

#define THICKNESS .05 

#define DIFF(X, T) DiffCh10((X), (T))

#define HEAT_MODEL
#define MASS_MODEL

matrix* CreateElementMatrix(struct fe1d *, Elem1D *, matrix *);
matrix* CreateDTimeMatrix(struct fe1d *, Elem1D *, matrix *);
matrix* CreateElementLoad(struct fe1d *, Elem1D *, matrix *);
int IsOnRightBoundary(struct fe1d *, int);
int IsOnLeftBoundary(struct fe1d *, int);
void ApplyAllBCs(struct fe1d *);
double DeformationGrad(struct fe1d *, double, double);

#endif

