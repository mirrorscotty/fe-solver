/**
 * @file integrate.c
 * Contains all sorts of stuff to integrate the interpolation functions.
 * Stuff that isn't properly documented is probably crap.
 */

#include <stdio.h>
#include <math.h>
#include "integrate.h"
#include "isoparam.h"
#include "basis.h"

/**
 * The number of Gauss points to use during integration. Change to 5 for
 * slightly more accurate answers when integrating weird, nonlinear functions.
 * You'll also need to change x3 and w3 to x5 and w5, respectively.
 */
#define NPTS 3 

/**
 * Coordinates of the Guass points when using 3 point Gaussian quadtrature
 */
double x3[] = {-0.774596669241483,
                0.,
                0.774596669241483};
/**
 * Weights for each point. (3 point)
 */
double w3[] = {0.555555555555555,
               0.888888888888888,
               0.555555555555555};

/**
 * Coordinates of the Gauss points (5 point)
 */
double x5[] = {-0.906179845938664,
               -0.538469310105683,
                0,
                0.538469310105683,
                0.906179845938664};
/**
 * Weights for each point (5 point)
 */
double w5[] = {0.236926885056189,
               0.478628670499366,
               0.568888888888889,
               0.478628670499366,
               0.236926885056189};


/* Integrate a function f from 0 to 1 using 3 point Gaussian Quadrature */
double quad3( double (*f)(double) )
{
    double result = 0;
    int i;

    double *x, *w;
    x = x5;
    w = w5;

    for(i=0; i<NPTS; i++) {
        result += w[i]/2 * (*f)( (x[i]+1)/2 );
    }

    return result;
}

/* Integrate a function f from 0 to 1 using 3 point Gaussian Quadrature */
double quad300(double (*f)(double), double (*g)(double))
{
    int i;
    double result = 0;
    double *x, *w;
    x = x5;
    w = w5;

    for(i=0; i<NPTS; i++) {
        result += w[i]/2 * (*f)( (x[i]+1)/2) * (*g)( (x[i]+1)/2);
    }

    return result;
}

/* Integrate a function f from 0 to 1 using 3 point Gaussian Quadrature */
double quad311(double (*f)(double), double (*g)(double))
{
    int i;
    double result = 0;

    double *x, *w;
    x = x5;
    w = w5;

    for(i=0; i<NPTS; i++) {
        result += w[i]/2 * diff(f, .5*(x[i]+1)) * diff(g, .5*(x[i]+1));
    }

    return result;
}

/* Integrate a function f from 0 to 1 using 3 point Gaussian Quadrature */
double quad301(double (*f)(double), double (*g)(double))
{
    int i;
    double result = 0;

    double *x, *w;
    x = x5;
    w = w5;

    for(i=0; i<NPTS; i++) {
        result += w[i]/2 * (*f)( .5*(x[i]+1)) * diff(g, .5*(x[i]+1));
    }

    return result;
}
/* Calculate the derivative of a function at point x using a centered
 * difference formula */
double diff(double (*f)(double), double x)
{
    double h = 1e-14;
    return (f(x+h)-f(x-h))/(2*h);
}

//double quad2d3(double (*f)(basis*,double,double))
//{
//    int i, j;
//    double result = 0;
//    double *x, *w;
//    x = x5;
//    w = w5;

//  for(i=0; i<NPTS; i++) {
//        for(j=0; j<NPTS; j++) {
//            result += w[i]/2 * w[j]/2 * (*f)( (x[i]+1)/2 , (x[j]+1)/2 );
// }
//    }
//    return result;
//}

/* Numerically integrate the desired basis function using Gaussian Quadrature.
 * This function only integrates 2D basis functions, and requires that the
 * element where the integration is being performed be passed as an argument
 * so that the result can be isoparametrically mapped correctly.
 */
double quad2d3(struct fe *p, Elem2D *elem, int func, int dx, int dy)
{
    int i, j;
    basis *b;
    double result = 0;
    double *x, *w;
    x = x3;
    w = w3;
    b = p->b;

    for(i=0; i<NPTS; i++) {
        for(j=0; j<NPTS; j++) {
            if(dx == 1) {
                result += w[i]/2 * w[j]/2
                    * ( b->Eval2Dx(b, func, (x[i]+1)/2, (x[j]+1)/2)
                    * IMapYEta(p, elem, (x[i]+1)/2, (x[j]+1)/2)
                    - b->Eval2Dy(b, func, (x[i]+1)/2, (x[j]+1)/2)
                    * IMapYXi(p, elem, (x[i]+1)/2, (x[j]+1)/2) )
                    * 1/IMapJ(p, elem, (x[i]+1)/2, (x[j]+1)/2)
                    * IMapCyl(p, elem, (x[i]+1)/2, (x[j]+1)/2);


            } else if(dy == 1) {
                result += w[i]/2 * w[j]/2 
                    * ( b->Eval2Dy(b, func, (x[i]+1)/2, (x[j]+1)/2)
                    * IMapXXi(p, elem, (x[i]+1)/2, (x[j]+1)/2)
                    - b->Eval2Dx(b, func, (x[i]+1)/2, (x[j]+1)/2)
                    * IMapXEta(p, elem, (x[i]+1)/2, (x[j]+1)/2) )
                    * 1/IMapJ(p, elem, (x[i]+1)/2, (x[j]+1)/2)
                    * IMapCyl(p, elem, (x[i]+1)/2, (x[j]+1)/2);
                
            } else {
                result += w[i]/2 * w[j]/2 
                    * b->Eval2D(b, func, (x[i]+1)/2, (x[j]+1)/2)
                    * 1/IMapJ(p, elem, (x[i]+1)/2, (x[j]+1)/2)
                    * IMapCyl(p, elem, (x[i]+1)/2, (x[j]+1)/2);
                
            }
            printf("i = %d, j = %d, result = %g, XXi = %g, XEta = %g, J = %g\n",
                    i, j, result,
                    IMapXXi(p, elem, (x[i]+1)/2, (x[j]+1)/2),
                    IMapXEta(p, elem, (x[i]+1)/2, (x[j]+1)/2),
                    IMapJ(p, elem, (x[i]+1)/2, (x[j]+1)/2));
        }
    }

    return result;
}

/* TODO: Make this function call more consise */
double quad2d32d3(struct fe *p, Elem2D *elem, int func1, int func2, int dx, int dy)
{
    int i, j;
    basis *b;
    double result = 0;
    double tmp1 = 0, tmp2 = 0;
    double *x, *w;
    x = x3;
    w = w3;
    b = p->b;

    for(i=0; i<NPTS; i++) {
        for(j=0; j<NPTS; j++) {
            if(dx == 1) {
                tmp1 += 1
                    * ( b->Eval2Dx(b, func1, (x[i]+1)/2, (x[j]+1)/2)
                    * IMapYEta(p, elem, (x[i]+1)/2, (x[j]+1)/2)
                    - b->Eval2Dy(b, func1, (x[i]+1)/2, (x[j]+1)/2)
                    * IMapYXi(p, elem, (x[i]+1)/2, (x[j]+1)/2) )
                    * 1/IMapJ(p, elem, (x[i]+1)/2, (x[j]+1)/2);
                tmp2 += 1
                    * ( b->Eval2Dx(b, func2, (x[i]+1)/2, (x[j]+1)/2)
                    * IMapYEta(p, elem, (x[i]+1)/2, (x[j]+1)/2)
                    - b->Eval2Dy(b, func2, (x[i]+1)/2, (x[j]+1)/2)
                    * IMapYXi(p, elem, (x[i]+1)/2, (x[j]+1)/2) )
                    * 1/IMapJ(p, elem, (x[i]+1)/2, (x[j]+1)/2);

            } else if(dy == 1) {
                tmp1 += 1
                    * ( b->Eval2Dy(b, func1, (x[i]+1)/2, (x[j]+1)/2)
                    * IMapXXi(p, elem, (x[i]+1)/2, (x[j]+1)/2)
                    - b->Eval2Dx(b, func1, (x[i]+1)/2, (x[j]+1)/2)
                    * IMapXEta(p, elem, (x[i]+1)/2, (x[j]+1)/2) )
                    * 1/IMapJ(p, elem, (x[i]+1)/2, (x[j]+1)/2);
                
                tmp2 += 1
                    * ( b->Eval2Dy(b, func2, (x[i]+1)/2, (x[j]+1)/2)
                    * IMapXXi(p, elem, (x[i]+1)/2, (x[j]+1)/2)
                    - b->Eval2Dx(b, func2, (x[i]+1)/2, (x[j]+1)/2)
                    * IMapXEta(p, elem, (x[i]+1)/2, (x[j]+1)/2) )
                    * 1/IMapJ(p, elem, (x[i]+1)/2, (x[j]+1)/2);
                
            } else {
                tmp1 += 1
                    * b->Eval2D(b, func1, (x[i]+1)/2, (x[j]+1)/2)
                    * 1/IMapJ(p, elem, (x[i]+1)/2, (x[j]+1)/2)
                    * IMapCyl(p, elem, (x[i]+1)/2, (x[j]+1)/2);
                
                tmp2 += 1
                    * b->Eval2D(b, func2, (x[i]+1)/2, (x[j]+1)/2)
                    * 1/IMapJ(p, elem, (x[i]+1)/2, (x[j]+1)/2)
                    * IMapCyl(p, elem, (x[i]+1)/2, (x[j]+1)/2);
            }
            result += w[i]/2 * w[j]/2 * tmp1*tmp2 * IMapCyl(p, elem, (x[i]+1)/2, (x[j]+1)/2);
            tmp1 = 0;
            tmp2 = 0;
            //printf("i = %d, j = %d, result = %g, XXi = %g, XEta = %g, J = %g\n",
            //        i, j, result,
            //        IMapXXi(p, elem, (x[i]+1)/2, (x[j]+1)/2),
            //        IMapXEta(p, elem, (x[i]+1)/2, (x[j]+1)/2),
            //        IMapJ(p, elem, (x[i]+1)/2, (x[j]+1)/2));
        }
    }

    return result;
}

/**
 * @brief Integrate a (1d) weighted residual so that the results can be stuck
 * into a local element matrix.
 * @param p The problem definition structure
 * @param guess A matrix containing the values of the residual at each of the
 * nodes. Row 1 of the matrix corresponds to node 1, etc. This is really only
 * used for the nonlinear solver. For the linear solver, this will be a matrix
 * of zeros.
 * @param residual A function pointer to the residual to integrate. The first
 * argument should be the same FE struct that was passed to this function. The
 * second argument is the guess matrix, the third is the pointer to the element
 * of interest, the fourth is the x coordinate (in isoparametric coordinates),
 * and the fifth and sixth are the row and column, respectively, in the element
 * matrix.
 * @param f1 The row of the element matrix
 * @param f2 The column of the element matrix
 * @return The value of the integrated function.
 */
double quad1d3generic(struct fe1d *p, matrix *guess, Elem1D *elem,
                      double (*residual)(struct fe1d*, matrix*, Elem1D*, double, int, int),
                      int f1, int f2)
{
    int i;
    double result = 0;
    double *x, *w;
    x = x3;
    w = w3;

    for(i=0; i<NPTS; i++) {
        result += w[i]/2 * residual(p, guess, elem, (x[i]+1)/2, f1, f2);
    }

    return result;
}

/**
 * The same as the 1D version, except the weighted residual now needs both the x
 * and y-coordinates.
 * @see quad1d3generic
 */
double quad2d3generic(struct fe *p, matrix *guess, Elem2D *elem,
                      double (*residual)(struct fe*, matrix*, Elem2D*, double, double, int, int),
                      int f1, int f2)
{
    int i, j;
    double result = 0;
    double *x, *w;
    x = x3;
    w = w3;

    for(i=0; i<NPTS; i++) {
        for(j=0; j<NPTS; j++) {
            result += w[i]/2 * w[j]/2 * residual(p, guess, elem, (x[i]+1)/2, (x[j]+1)/2, f1, f2);
        }
    }

    return result;
}

double quad2d3arc(struct fe *p, matrix *guess, matrix *prevguess, Elem2D *elem,
                  double (*residual)(struct fe*, matrix*, matrix*, Elem2D*, double, double))
{
    int i, j;
    double result = 0;
    double *x, *w;
    x = x3;
    w = w3;

    for(i=0; i<NPTS; i++) {
        for(j=0; j<NPTS; j++) {
            result += w[i]/2 * w[j]/2 * residual(p, guess, prevguess, elem, (x[i]+1)/2, (x[j]+1)/2);
        }
    }

    return result;
}

double quad2d3tri(struct fe *p, matrix *guess, Elem2D *elem,
                      double (*residual)(struct fe*, matrix*, Elem2D*, double, double, int, int),
                      int f1, int f2)
{
    int i, j;
    double result = 0;
    double x[] = { 0.35355339, 0.70710678 };
    double y[] = { 0.35355339, 0.70710678 };
    double w[] = { -9/32, 25/96, 25/96, 25/96 };

    for(i=0; i<2; i++) {
        for(j=0; j<2; j++) {
            result += w[i*j] * residual(p, guess, elem, i, j, x[i], y[j]);
        }
    }
    return result;
}

double diff2dx(basis *b, int func, double x, double y)
{
    double h = 1e-14;
    return (b->Eval2D(b, func, x+h, y) - b->Eval2D(b, func, x-h, y))/(2*h);
}

double diff2dy(basis *b, int func, double x, double y)
{
    double h = 1e-14;
    return (b->Eval2D(b, func, x, y+h) - b->Eval2D(b, func, x, y-h))/(2*h);
}
