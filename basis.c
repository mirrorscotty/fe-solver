#include <math.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include "basis.h"
#include "integrate.h"

/* //Test function
int main(int argc, char *argv[])
{
    basis *b;
    b = MakeLinBasis(2);

    printf("p1 = %g\np2 = %g\np3 = %g\np4 = %g\n",
           quad2d3(b, 0, 0, 1),
           quad2d3(b, 1, 0, 1),
           quad2d3(b, 2, 0, 1),
           quad2d3(b, 3, 0, 1));
          // EvalLin2D(b, 0, 0, 0),
          // EvalLin2D(b, 1, 0, 0),
          // EvalLin2D(b, 2, 0, 0),
          // EvalLin2D(b, 3, 0, 0));

    DestroyBasis(b);
    return 42;
}
*/

/* Make a set of linear basis functions. Includes information on how to
   properly assemble the global matricies. */
basis* MakeLinBasis(int dimension)
{
    basis *b;
    b = NULL;
    b = (basis*) malloc(sizeof(basis));
    if(!b) {
        return NULL;
    }
    b->phi = b->dphi = NULL;

    b->n = 2 * dimension;
    b->overlap = 1 * dimension;

    b->phi = (double(**)(double)) malloc(sizeof(double(*)(double))*b->n);
    b->dphi = (double(**)(double)) malloc(sizeof(double(*)(double))*b->n);
    if(b->phi) {
        b->phi[0] = &lin1d1;
        b->phi[1] = &lin1d2;
    } else {
        fprintf(stderr, "Failed to allocate %lu bytes.\n", sizeof(double(*)(double))*b->n);
    }

    if(b->dphi) {
        b->dphi[0] = &dlin1d1;
        b->dphi[1] = &dlin1d2;
    } else {
        fprintf(stderr, "Failed to allocate %lu bytes.\n", sizeof(double(*)(double))*b->n);
    }

    b->dim = dimension;

    return b;
}

basis* MakeQuadBasis(int dimension)
{
    basis *b;
    b = NULL;
    b = (basis*) malloc(sizeof(basis));
    if(!b) {
        return NULL;
    }
    b->phi = b->dphi = NULL;

    b->n = 3 * dimension;
    b->overlap = 1 * dimension;

    b->phi = (double(**)(double)) calloc(b->n, sizeof(double(**)(double)));
    b->dphi = (double(**)(double)) calloc(b->n, sizeof(double(**)(double)));
    if(b->phi) {
        b->phi[0] = &quad1d1;
        b->phi[1] = &quad1d2;
        b->phi[2] = &quad1d3;
    } else {
        fprintf(stderr, "Failed to allocate %lu bytes.\n", sizeof(double(*)(double))*b->n);
    }

    if(b->dphi) {
        b->dphi[0] = &dquad1d1;
        b->dphi[1] = &dquad1d2;
        b->dphi[2] = &dquad1d3;
    } else {
        fprintf(stderr, "Failed to allocate $lu bytes.\n", sizeof(double(*)(double))*b->n);
    }

    b->dim = dimension;

    return b;
}

basis* MakeCubicBasis(int dimension)
{
    basis *b;
    b = NULL;
    b = (basis*) malloc(sizeof(basis));
    if(!b) {
        return NULL;
    }
    b->phi = b->dphi = NULL;

    b->n = 4 * dimension;
    b->overlap = 2 * dimension;

    b->phi = (double(**)(double)) malloc(sizeof(double(*)(double))*b->n);
    b->dphi = (double(**)(double)) malloc(sizeof(double(*)(double))*b->n);
    if(b->phi) {
        b->phi[0] = &cubic1d1;
        b->phi[1] = &cubic1d2;
        b->phi[2] = &cubic1d3;
        b->phi[3] = &cubic1d4;
    } else {
        fprintf(stderr, "Failed to allocate %lu bytes.\n", sizeof(double(*)(double))*b->n);
    }

    if(b->dphi) {
        b->dphi[0] = &dcubic1d1;
        b->dphi[1] = &dcubic1d2;
        b->dphi[2] = &dcubic1d3;
        b->dphi[3] = &dcubic1d4;
    } else {
        fprintf(stderr, "Failed to allocate %lu bytes.\n", sizeof(double(*)(double))*b->n);
    }

    b->dim = dimension;

    return b;
}
    
/* Delete stuff */
void DestroyBasis(basis *b)
{
    if(b->phi)
        free(b->phi);
    if(b->dphi)
        free(b->dphi);
    if(b)
        free(b);
    return;
}

/* Evaluate a linear 2 dimensional basis function. The second argument
   determines which of the bilinear basis functions to evaluate. The last two
   arguments are the x and y coordinates in isoparametric space to evaluate the
   function at. */
double EvalLin2D(basis *b, int func, double x, double y)
{
    //printf("Evaluating function (%d, %d)\n", (func>1)?1:0, func%2);
    int i, j;
    /* Determine which of the linear basis functions to use */
    i = (func>1)?1:0;
    j = func%2;
    return b->phi[i](x) * b->phi[j](y);
    //return EvalBasis(b, (func>1)?1:0, x, func % 2, y);
}

double EvalLin2Dx(basis *b, int func, double x, double y)
{
    int i, j;
    i = (func>1)?1:0;
    j = func%2;
    return b->dphi[i](x) * b->phi[j](y);
}

double EvalLin2Dy(basis *b, int func, double x, double y)
{
    int i, j;
    i = (func>1)?1:0;
    j = func%2;
    return b->phi[i](x) * b->dphi[j](y);
}

/* Not even slightly finished yet. */
//double EvalBasisV(basis *b, int function, vector *v)
//{
//    int i;

//    if(len(v) != b->dim) {
//        fprintf(stderr, "Cataclysmic failure! Exiting immediately!\n");
//        exit(0);
//    }

//    for(i=0; i<1;i++)
//}

/* General function for evaluating a basis at a specified point. The first
   argument is the basis to evaluate, and the next two arguments are the
   1d function number and the isoparametric coordinate in that dimension. For
   example, EvalBasis(b, 0, .5, 0, .2) would evaluate the first 2d bilinear
   basis function at (.5, .2). The dimension of the basis must be equal to the
   number of coordinates supplied. */
double EvalBasis(basis *b, ... )
{
    va_list arglist;
    int *phi;
    double *arg;
    int i;
    double result = 1;

    va_start(arglist, 2*b->dim);
    phi = (int*) calloc(b->dim, sizeof(int));
    arg = (double*) calloc(b->dim, sizeof(double));

    for(i=0; i<b->dim; i++) {
        phi[i] = va_arg(arglist, int);
        //printf("%d\n", phi[i]);
    }

    for(i=0; i<b->dim; i++) {
        arg[i] = va_arg(arglist, double);
        //printf("%g\n", arg[i]);
    }

    for(i=0; i<b->dim; i++) {
        result *= b->phi[i](arg[i]);
        //printf("%g\n", result);
    }

    va_end(arglist);
    free(arg);
    free(phi);
    return result;
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
