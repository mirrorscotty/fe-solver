/**
 * @file heat-exact.c
 * Exact solution to the heat equation. For use as a comparison against the FEM
 * solution to validate the numerical model.
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "matrix.h"

struct heat {
    double L; /* Domain width */
    double (*f)(double); /* Initial Condition */
    double a; /* Thermal diffusivity */
    int nmax; /* Number of terms of the solution to use */
};

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

double Init(double x)
{
    return sin(M_PI*x);
}

double InitQuad(double x)
{
    return -(x-1)*(x-1)+1;
}

/* Don't use */
double Dn_gauss(struct heat *h, int n)
{
    int i;
    double result = 0;
    double L = h->L;

    for(i=0; i<5; i++)
        result += w5[i] * h->f(L*(x5[i]+1)/2) * sin(n*M_PI/2*(x5[i]+1));
        //result += w3[i] *  h->f(h->L/2*x3[i]+h->L/2)
        //          * sin(n*M_PI/h->L * (h->L/2*x3[i]+h->L/2));

    //result *= 2/h->L;

    return result;
}

/* Also broken */
double Dn(struct heat *h, int n)
{
    double result = 0;
    double L = h->L;
    double pi = M_PI;

    printf("L = %g\npi = %g\nn = %d\n", L, pi, n);

    return (2*n*sin(pi*L)*cos(n*pi)-2*L*cos(pi*L)*sin(pi*n))/(pi*(L-n)*(L+n));
}

/* Uses the initial condition f(x) = -(x-1)^2 + 1
 * To use this function, L must be 2.
 */
double Dn_quad(struct heat *h, int n)
{
    double result = 0;
    double L = h->L;
    double pi = M_PI;

    return (L*L*((pi*pi*(L-2)*n*n-2*L)*cos(pi*n)-2*pi*(L-1)*n*sin(pi*n)+2*L))/(pi*pi*pi*n*n*n);
}

double heat_exact(struct heat *h, double x, double t)
{
    int i;
    double result = 0;

    for(i=1; i<=h->nmax; i++)
    {
        result += Dn_quad(h, i) * sin(i*M_PI*x/h->L)
            * exp( -i*i * M_PI*M_PI * h->a * t / (h->L*h->L) );
    }

    return result;
}

/* Initialize a "heat" structure and set the default values for the problem. */
struct heat* heat_default()
{
    struct heat *h;
    h = (struct heat *) calloc(1, sizeof(struct heat));

    h->a = 1;
    h->L = 2;
    h->f = &InitQuad; /* Doesn't really mattter. */
    h->nmax = 100;

    return h;
}

int main(int argc, char *argv[])
{
    int i;

    double t = 3; /* The time to output the solution at. */
    int npts = 6; /* The number of points to output. */
    struct heat *h;
    matrix *answer;

    h = heat_default();
    answer = CreateMatrix(npts, 1);

    for(i=0; i<npts; i++) {
        setval(answer, heat_exact(h, h->L/(npts-1)*i, t), i, 0);
    }

    mtxprnt(answer);
    free(h);
    DestroyMatrix(answer);

    return 0;
}

