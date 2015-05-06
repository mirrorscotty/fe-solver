#ifndef COMMON_H
#define COMMON_H

//#define TVAR 0
#define CVAR 0 
#define SVAR 1
//#define PVAR 1

//#define TAMB 310 // K
//#define TAMB 500 // K
//#define TINIT 273 //K
//#define TINIT 313.15 //K
//#define TINIT 353.15 //K
#define TINIT 333 //K
#define HCONV 50

#define CAMB 0.15 // kg/kg db
//#define CAMB 0.0861867 // kg/kg db (RH=.65, T=80C)
//#define CAMB 0.138 // kg/kg db (RH=.65, T=40C)
//#define CAMB 0.0975028 // kg/kg db (RH=.7, T=80C)
//#define CAMB 0.15049 // kg/kg db (RH=.7, T=40C)
#define CINIT .5 // kg/kg db
#define KC_CONV 2e-11

#define THICKNESS 1e-3

#define POISSON .37

//#define SHRINKAGE

#define DIFF(X, T) DiffCh10((X), (T))
//#define DIFF(X, T) DiffAvg((X), (T))

#define CREEP(TIME, T, X, P, DERIV)  CreepLaura2((TIME), (T), (X), (P), (DERIV))

#define BETA .6

double DiffAvg(double, double);

matrix* CreateElementMatrix(struct fe1d *, Elem1D *, matrix *);
matrix* CreateDTimeMatrix(struct fe1d *, Elem1D *, matrix *);
matrix* CreateElementLoad(struct fe1d *, Elem1D *, matrix *);
int IsOnRightBoundary(struct fe1d *, int);
int IsOnLeftBoundary(struct fe1d *, int);
void ApplyAllBCs(struct fe1d *);
//double DeformationGrad(struct fe1d *, double, double);

#endif

