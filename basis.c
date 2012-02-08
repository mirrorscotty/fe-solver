#include <math.h>
#include <stdlib.h>
#include "basis.h"
#include "integrate.h"

basis* MakeLinBasis()
{
    basis *b;
    b = (basis*) malloc(sizeof(basis*));

    b->n = 2;
    b->overlap = 1;

    b->phi = (double(**)(double)) malloc(sizeof(double(*)(double))*b->n);
    b->dphi = (double(**)(double)) malloc(sizeof(double(*)(double))*b->n);
    b->phi[0] = &lin1d1;
    b->phi[1] = &lin1d2;

    b->dphi[0] = &dlin1d1;
    b->dphi[1] = &dlin1d2;

    return b;
}

basis* MakeQuadBasis()
{
    basis *b;
    b = (basis*) malloc(sizeof(basis*));

    b->n = 3;
    b->overlap = 1;

    b->phi = (double(**)(double)) malloc(sizeof(double(*)(double))*b->n);
    b->dphi = (double(**)(double)) malloc(sizeof(double(*)(double))*b->n);
    b->phi[0] = &quad1d1;
    b->phi[1] = &quad1d2;
    b->phi[2] = &quad1d3;

    b->dphi[0] = &dquad1d1;
    b->dphi[1] = &dquad1d2;
    b->dphi[2] = &dquad1d3;

    return b;
}

basis* MakeCubicBasis()
{
    basis *b;
    b = (basis*) malloc(sizeof(basis*));

    b->n = 4;
    b->overlap = 2;

    b->phi = (double(**)(double)) malloc(sizeof(double(*)(double))*b->n);
    b->dphi = (double(**)(double)) malloc(sizeof(double(*)(double))*b->n);
    b->phi[0] = &cubic1d1;
    b->phi[1] = &cubic1d2;
    b->phi[2] = &cubic1d3;
    b->phi[3] = &cubic1d4;

    b->dphi[0] = &dcubic1d1;
    b->dphi[1] = &dcubic1d2;
    b->dphi[2] = &dcubic1d3;
    b->dphi[3] = &dcubic1d4;

    return b;
}
    
void DestroyBasis(basis *b)
{
    free(b->phi);
    free(b->dphi);
    free(b);
    return;
}

/* 1d Linear Basis Functions */
double lin1d1(double x) { return 1-x; }
double lin1d2(double x) { return x; }

double dlin1d1(double x) { return -1; }
double dlin1d2(double x) { return 1; }

/* 1d Quadratic Basis Functions */
double quad1d1(double x) { return 1-3*x+2*pow(x, 2); }
double quad1d2(double x) { return 4*(x-pow(x, 2)); }
double quad1d3(double x) { return -x+2*pow(x, 2); }

double dquad1d1(double x) { return 4*x-3; }
double dquad1d2(double x) { return 4-8*x; }
double dquad1d3(double x) { return 4*x-1; }

/* 1d Hermite Cubic Functions */
double cubic1d1(double x) { return 1-3*pow(x,2)+2*pow(x,3); }
double cubic1d2(double x) { return x-2*pow(x,2)+pow(x,3); } 
double cubic1d3(double x) { return 3*pow(x,2)-2*pow(x,3); } 
double cubic1d4(double x) { return pow(x,3)-pow(x,2); }

double dcubic1d1(double x) { return 6*x*x-6*x; }
double dcubic1d2(double x) { return 1-4*x+3*x*x; }
double dcubic1d3(double x) { return 6*x-6*x*x; }
double dcubic1d4(double x) { return 3*x*x-2*x; }
