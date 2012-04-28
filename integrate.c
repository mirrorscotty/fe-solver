#include <stdio.h>
#include <math.h>
#include "integrate.h"
#include "isoparam.h"
#include "basis.h"

double square(double x) {
	return x*x*x;
}
//void main() {
//	printf("%1.10f", quad300(&square, &square));
//}
//
#define NPTS 5 // Number of Gauss points to use

double x3[] = {-0.774596669241483,
                0.,
                0.774596669241483};
double w3[] = {0.555555555555555,
               0.888888888888888,
               0.555555555555555};

double x5[] = {-0.906179845938664,
               -0.538469310105683,
                0,
                0.538469310105683,
                0.906179845938664};
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
    x = x5;
    w = w5;
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
