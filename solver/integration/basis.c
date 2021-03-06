/**
 * @file basis.c
 * Contains all the functions and structures for storing the basis functions
 * used in the finite element interpolation. Currently has stuff for linear,
 * quadratic, and hermite cubic functions. It also does 1D and 2D
 * interpolation.
 */

#include <math.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include "basis.h"
#include "integrate.h"

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
/**
 * @brief Make a set of linear basis functions.
 * Make a set of linear basis functions. Includes information on how to
 * properly assemble the global matricies.
 * @param dimension Dimension of the domain
 * @return Data structure with pointers to the basis functions
 */
basis* MakeLinBasis(int dimension)
{
    basis *b;
    b = NULL;
    b = (basis*) malloc(sizeof(basis));
    if(!b) {
        return NULL;
    }
    b->phi = b->dphi = NULL;

    b->n = pow(2, dimension);
    b->overlap = pow(1, dimension);

    b->phi = (double(**)(double)) malloc(sizeof(double(*)(double))*b->n);
    b->dphi = (double(**)(double)) malloc(sizeof(double(*)(double))*b->n);
    if(b->phi) {
        b->phi[0] = &lin1d1;
        b->phi[1] = &lin1d2;
    } else {
        fprintf(stderr, "Failed to allocate %lu bytes.\n", (unsigned long) sizeof(double(*)(double))*b->n);
    }

    if(b->dphi) {
        b->dphi[0] = &dlin1d1;
        b->dphi[1] = &dlin1d2;
    } else {
        fprintf(stderr, "Failed to allocate %lu bytes.\n", (unsigned long) sizeof(double(*)(double))*b->n);
    }

    b->dim = dimension;

    /* Make sure the basis can evaluate itself properly */
    if(b->dim == 2) {
        b->Eval2D = &EvalLin2D;
        b->Eval2Dx = &EvalLin2Dx;
        b->Eval2Dy = &EvalLin2Dy;
    } else {
        b->Eval2D = NULL;
        b->Eval2Dx = NULL;
        b->Eval2Dy = NULL;
    }

    return b;
}

/**
 * @brief Make a set of quadratic basis functions.
 * Make a set of quadratic basis functions. Includes information on how to
 * properly assemble the global matricies.
 * @param dimension Dimension of the domain
 * @return Data structure with pointers to the basis functions
 */
basis* MakeQuadBasis(int dimension)
{
    basis *b;
    b = NULL;
    b = (basis*) malloc(sizeof(basis));
    if(!b) {
        return NULL;
    }
    b->phi = b->dphi = NULL;

    b->n = pow(3, dimension);
    b->overlap = pow(1, dimension);

    b->phi = (double(**)(double)) calloc(b->n, sizeof(double(**)(double)));
    b->dphi = (double(**)(double)) calloc(b->n, sizeof(double(**)(double)));
    if(b->phi) {
        b->phi[0] = &quad1d1;
        b->phi[1] = &quad1d2;
        b->phi[2] = &quad1d3;
    } else {
        fprintf(stderr, "Failed to allocate %lu bytes.\n", (unsigned long) sizeof(double(*)(double))*b->n);
    }

    if(b->dphi) {
        b->dphi[0] = &dquad1d1;
        b->dphi[1] = &dquad1d2;
        b->dphi[2] = &dquad1d3;
    } else {
        fprintf(stderr, "Failed to allocate %lu bytes.\n", (unsigned long) sizeof(double(*)(double))*b->n);
    }

    b->dim = dimension;

    if(b->dim == 2) {
        b->Eval2D = &EvalQuad2D;
        b->Eval2Dx = &EvalQuad2Dx;
        b->Eval2Dy = &EvalQuad2Dy;
    } else {
        b->Eval2D = NULL;
        b->Eval2Dx = NULL;
        b->Eval2Dy = NULL;
    }

    return b;
}

/**
 * @brief Make a set of hermite cubic basis functions.
 * Make a set of cubic basis functions. Includes information on how to
 * properly assemble the global matricies.
 * @param dimension Dimension of the domain
 * @return Data structure with pointers to the basis functions
 */
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
        fprintf(stderr, "Failed to allocate %lu bytes.\n", (unsigned long) sizeof(double(*)(double))*b->n);
    }

    if(b->dphi) {
        b->dphi[0] = &dcubic1d1;
        b->dphi[1] = &dcubic1d2;
        b->dphi[2] = &dcubic1d3;
        b->dphi[3] = &dcubic1d4;
    } else {
        fprintf(stderr, "Failed to allocate %lu bytes.\n", (unsigned long) sizeof(double(*)(double))*b->n);
    }

    b->dim = dimension;

    return b;
}

basis* Make2DTriBasis()
{
    basis *b;
    b = NULL;
    b = (basis*) calloc(1, sizeof(basis));

    if(!b)
        return NULL;

    b->phi = b->dphi = NULL;
    b->n = 3;
    b->overlap = 42;

    b->phi2d = (double (**)(double, double)) calloc(b->n, sizeof(double(*)(double, double)));

    b->phi2d[1] = &lintri2d1;
    b->phi2d[2] = &lintri2d2;
    b->phi2d[3] = &lintri2d3;

    b->dim = 2;


    b->Eval2D = &EvalLinTri2D;
    b->Eval2Dx = &EvalLinTri2Dx;
    b->Eval2Dy = &EvalLinTri2Dy;

    return b;
}


/**
 * @brief Delete a basis!
 * @param b The basis to dispose of.i
 */
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

/**
 * @brief Evaluate a linear 2 dimensional basis function.
 *
 * The second argument
 * determines which of the bilinear basis functions to evaluate. The last two
 * arguments are the x and y coordinates in isoparametric space to evaluate the
 * function at.
 * @param b The (2d) basis containing the functions to evaluate
 * @param func Which function to evaluate
 * @param x The x-coordinate inside the element
 * @param y The y-coordinate
 * @return The result
 */
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

/**
 * @brief Evaluates the derivative of the basis functions with respect to x.
 * @param b The basis struct
 * @param func The function to evaluate
 * @param x X-coordinate
 * @param y Y-coordinate
 * @return The value at those coordinates
 */
double EvalLin2Dx(basis *b, int func, double x, double y)
{
    int i, j;
    i = (func>1)?1:0;
    j = func%2;
    return b->dphi[i](x) * b->phi[j](y);
}

/**
 * @brief Evaluates the derivative of the basis functions with respect to x.
 * @param b The basis struct
 * @param func The function to evaluate
 * @param x X-coordinate
 * @param y Y-coordinate
 * @return The value at those coordinates
 */
double EvalLin2Dy(basis *b, int func, double x, double y)
{
    int i, j;
    i = (func>1)?1:0;
    j = func%2;
    return b->phi[i](x) * b->dphi[j](y);
}

/**
 * @brief Evaluate a biquadritic basis
 * @see EvalLin2D
 */
double EvalQuad2D(basis *b, int func, double x, double y)
{
    int i, j;
    i = (func>5)?2:((func>2)?1:0);
    j = func%3;
    return b->phi[i](x) * b->phi[j](y);
}

/**
 * @brief Evaluate the x derivative of a set of biquadratic basis functions.
 * @see EvalLin2Dx
 */
double EvalQuad2Dx(basis *b, int func, double x, double y)
{
    int i, j;
    i = (func>5)?2:((func>2)?1:0);
    j = func%3;
    return b->dphi[i](x) * b->phi[j](y);
}

/* @brief Evaluate the y derivative of a set of biquadratic basis functions.
 * @see EvalLin2Dy
 */
double EvalQuad2Dy(basis *b, int func, double x, double y)
{
    int i, j;
    i = (func>5)?2:((func>2)?1:0);
    j = func%3;
    return b->phi[i](x) * b->dphi[j](y);
}

double EvalLinTri2D(basis *b, int func, double x, double y)
{
    return b->phi2d[func](x, y);
}

double EvalLinTri2Dx(basis *b, int func, double x, double y)
{
    double h = 1e-14;
    return (b->phi2d[func](x+h, y) - b->phi2d[func](x-h, y))/(2*h);
}

double EvalLinTri2Dy(basis *b, int func, double x, double y)
{
    double h = 1e-14;
    return (b->phi2d[func](x, y+h) - b->phi2d[func](x, y-h))/(2*h);
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
/* Pretty sure this doesn't work.
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
*/

/* 1d Linear Basis Functions */
double lin1d1(double x) { return 1-x; }
double lin1d2(double x) { return x; }
/* First derivatives */
double dlin1d1(double x) { return -1; }
double dlin1d2(double x) { return 1; }

/* 1d Quadratic Basis Functions */
double quad1d1(double x) { return 1-3*x+2*pow(x, 2); }
double quad1d2(double x) { return 4*(x-pow(x, 2)); }
double quad1d3(double x) { return -x+2*pow(x, 2); }
/* First derivatives */
double dquad1d1(double x) { return 4*x-3; }
double dquad1d2(double x) { return 4-8*x; }
double dquad1d3(double x) { return 4*x-1; }

/* 1d Hermite Cubic Functions */
double cubic1d1(double x) { return 1-3*pow(x,2)+2*pow(x,3); }
double cubic1d2(double x) { return x-2*pow(x,2)+pow(x,3); }
double cubic1d3(double x) { return 3*pow(x,2)-2*pow(x,3); }
double cubic1d4(double x) { return pow(x,3)-pow(x,2); }
/* First derivatives */
double dcubic1d1(double x) { return 6*x*x-6*x; }
double dcubic1d2(double x) { return 1-4*x+3*x*x; }
double dcubic1d3(double x) { return 6*x-6*x*x; }
double dcubic1d4(double x) { return 3*x*x-2*x; }

/* 2d Linear Triangle Basis Functions (Isoparametric Coordinates) */
double lintri2d1(double x, double y) { return x; }
double lintri2d2(double x, double y) { return y; }
double lintri2d3(double x, double y) { return 1-x-y; }

/* 2d Quadratic Triangle Basis Functions (Isoparametric Coordinates) */
double quadtri2d1(double x, double y) { return x*(2*x-1); }
double quadtri2d2(double x, double y) { return y*(2*y-1); }
double quadtri2d3(double x, double y) { return (1-x-y)*(2*(1-x-y)-1); }
double quadtri2d4(double x, double y) { return 4*x*y; }
double quadtri2d5(double x, double y) { return 4*y*(1-x-y); }
double quadtri2d6(double x, double y) { return 4*(1-x-y)*x; }

