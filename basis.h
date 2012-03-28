#ifndef BASIS_H
#define BASIS_H

typedef struct {
    int n; /* Number of functions in the basis */
    int overlap; /* Amount of overlap when assembling the global matricies */
    double (**phi)(double); /* Array of basis functions */
    double (**dphi)(double); /* First derivatives of the basis functions */
    int dim; /* Dimension of the problem */
} basis;

basis* MakeLinBasis(int);
basis* MakeQuadBasis(int);
basis* MakeCubicBasis(int);

void DestroyBasis(basis*);

double EvalBasis(basis*, ... );

double lin1d1(double);
double lin1d2(double);
double dlin1d1(double);
double dlin1d2(double);

double quad1d1(double);
double quad1d2(double);
double quad1d3(double);
double dquad1d1(double);
double dquad1d2(double);
double dquad1d3(double);

double cubic1d1(double);
double cubic1d2(double);
double cubic1d3(double);
double cubic1d4(double);
double dcubic1d1(double);
double dcubic1d2(double);
double dcubic1d3(double);
double dcubic1d4(double);

#endif
