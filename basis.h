#ifndef BASIS_H
#define BASIS_H

struct _basis {
    int n; /* Number of functions in the basis */
    int overlap; /* Amount of overlap when assembling the global matricies */
    double (**phi)(double); /* Array of basis functions */
    double (**dphi)(double); /* First derivatives of the basis functions */
    double (**phi2d)(double, double); /* Only for 2d triangles */
    int dim; /* Dimension of the problem */

    double (*Eval2D)(struct _basis*, int, double, double);
    double (*Eval2Dx)(struct _basis*, int, double, double);
    double (*Eval2Dy)(struct _basis*, int, double, double);
};
typedef struct _basis basis;

basis* MakeLinBasis(int);
basis* MakeQuadBasis(int);
basis* MakeCubicBasis(int);

void DestroyBasis(basis*);

double EvalLin2D(basis*, int, double, double);
double EvalLin2Dx(basis*, int, double, double);
double EvalLin2Dy(basis*, int, double, double);
double EvalBasis(basis*, ... );

double EvalQuad2D(basis*, int, double, double);
double EvalQuad2Dx(basis*, int, double, double);
double EvalQuad2Dy(basis*, int, double, double);

double EvalLinTri2D(basis*, int, double, double);
double EvalLinTri2Dx(basis*, int, double, double);
double EvalLinTri2Dy(basis*, int, double, double);

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

double lintri2d1(double, double);
double lintri2d2(double, double);
double lintri2d3(double, double);

#endif
