/**
 * @file heat-gui.c
 * Set of functions to solve the following differential equations:
 * \f[
 * \frac{\partial T}{\partial t} = \alpha(T)\frac{\partial^2 T}{\partial t^2}
 * \f]
 * \f[
 * c_1 = A_1 \exp{-\frac{E_{a,1}}{R T}}
 * \f]
 * \f[
 * c_2 = A_2 \exp{-\frac{E_{a,2}}{R T}}
 * \f]
 * where \f$\alpha(T)\f$ is the temperature-dependent thermal diffusivity, and
 * \f$c_1\f$ and \f$c_2\f$ are the concentrations of two components that whose
 * rate of reaction depends on temperature.
 *
 * The heat conduction equation in this file is not correct. The appropriate
 * nonlinear heat conduction equation is:
 * \f[
 * \rho C_p\frac{\partial T}{\partial t} = \nabla\cdot\left(k\nabla T\right)
 * \f]
 * However, the form presented here should be somewhat accurate.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "matrix.h"
#include "basis.h"
#include "integrate.h"
#include "mesh1d.h"
#include "isoparam.h"
#include "finite-element1d.h"

#include "material-data/choi-okos/choi-okos.h"

#include "heat-transfer.h"

extern choi_okos *comp_global;

/* The following two functions are for heat conduction. */
/* Creates the Jacobian and helps solve for the current time step */
double Residual(struct fe1d *p, matrix *guess, Elem1D *elem, double x, int f1, int f2)
{
    double term1 = 0, term2 = 0;
    double h = 1e-7;
    basis *b;
    b = p->b;

    /* This is just to fetch the value of T */
    double T;
    solution *s;
    //s = FetchSolution(p, p->t-1);
    s = CreateSolution(p->t, p->dt, guess);
    //if(!s)
    if(p->t==0)
        T = scaleTemp(p->charvals, 273);
    else
        T = EvalSoln1D(p, 0, elem, s, x);
    /* Now that we've calculated T, we no longer need this */
    free(s);

    term1  = b->dphi[f1](x) * b->dphi[f2](x);
    term1 *= IMap1D(p, elem, x);
    //term1 *= (alpha(comp_global, uscaleTemp(p->charvals, T))/p->charvals.alpha)
    //            * b->phi[f1](x) * IMap1D(p, elem, x);
    term1 *= k(comp_global, uscaleTemp(p->charvals, T))/p->charvals.k;

    term2 = (k(comp_global, T+h)-k(comp_global, T-h))/(2*h);
    term2 *= 1/p->charvals.k;
    term2 *= pow(b->dphi[f1](x), 2) * b->phi[f2](x);
    term2 *= pow(IMap1D(p, elem, x), 2);

    return term1+term2;
}

/* Calculate the coefficient matrix for the time derivative unknowns */
double ResDt(struct fe1d *p, matrix *guess, Elem1D *elem, double x, int f1, int f2)
{
    double T;
    solution *s;

    double value;
    basis *b;
    b = p->b;

    s = CreateSolution(p->t, p->dt, guess);
    T = EvalSoln1D(p, 0, elem, s, x);
    /* Now that we've calculated T, we no longer need this */
    free(s);

    value = b->phi[f1](x) * b->phi[f2](x) * 1/IMap1D(p, elem, x);
    value *= rho(comp_global, uscaleTemp(p->charvals, T))
            * Cp(comp_global, uscaleTemp(p->charvals, T))
            * p->charvals.alpha/p->charvals.k;
    value *= IMapCyl1D(p, elem, x);

    return value;
}

matrix* CreateElementMatrix(struct fe1d *p, Elem1D *elem, matrix *guess)
{
    basis *b;
    b = p->b;
    
    int v = p->nvars;
    
    int i, j;
    double value = 0;
    matrix *m;
    
    m = CreateMatrix(b->n*v, b->n*v);
    
    for(i=0; i<b->n*v; i+=v) {
        for(j=0; j<b->n*v; j+=v) {
            value = quad1d3generic(p, guess, elem, &Residual, i/v, j/v);
            setval(m, value, i, j);
        }
    }

    return m;
}

/* Create the coefficient matrix for the time derivatives */
matrix* CreateDTimeMatrix(struct fe1d *p, Elem1D *elem, matrix *guess) {
    basis *b;
    b = p->b;
    
    int v = p->nvars;
    
    int i, j;
    double value = 0;
    matrix *m;
    
    m = CreateMatrix(b->n*v, b->n*v);
    
    for(i=0; i<b->n*v; i+=v) {
        for(j=0; j<b->n*v; j+=v) {
            value = quad1d3generic(p, guess, elem, &ResDt, i/v, j/v);
            setval(m, value, i, j);
        }
    }

    return m;
}

/* Create the load vector... thing */
matrix* CreateElementLoad(struct fe1d *p, Elem1D *elem, matrix *guess) {
    basis *b;
    b = p->b;
    
    int v = p->nvars;
    
    matrix *m;
    
    m = CreateMatrix(b->n*v, 1);

    return m;
}

/* Returns true if the specified row (node) is on the right-most boundary of the
 * domain. */
int IsOnRightBoundary(struct fe1d *p, int row)
{
    double width = p->mesh->x2 - p->mesh->x1;
    double x = valV(p->mesh->nodes, row/p->nvars);
    
    if(fabs(x - width) < 1e-5)
        return 1;
    else
        return 0;
}

/* Same as above, only for the left boundary */
int IsOnLeftBoundary(struct fe1d *p, int row)
{
    double x = valV(p->mesh->nodes, row/p->nvars);
  
    if(fabs(x) < 1e-5)
        return 1;
    else
        return 0;
}

double ExternalTemp(struct fe1d *p, int row)
{
    return scaleTemp(p->charvals, p->charvals.Tc);
}

/* The way this function is implemented is probably not mathematically accurate.
 * Strictly speaking, for the implicit solver, the temperature fetched should be
 * the temperature at the next time step, not the one that we already have the
 * solution for. The way it is now should be good enough (tm).*/
double ConvBC(struct fe1d *p, int row)
{
    double Tinf = ExternalTemp(p, row);
    double T;
    double Bi = BiotNumber(p->charvals);

    /* Fetch the value of T from the previous solution. If this is being
     * applied at the initial condition, then simply return 0.*/
    //solution *s;
    //s = FetchSolution(p, p->t-1);
    //if(s) {
        //T = val(s->val, row, 0);
    if(p->guess) {
        T = val(p->guess, row, 0);
        if(T==0)
            return 0;
        else
            return -Bi*(T-Tinf);
    } else {
        return 0;
    }
}

void ApplyAllBCs(struct fe1d *p)
{
    double Bi = BiotNumber(p->charvals); 
    
    /* BC at x=L
     * This approximates any Biot number larger than 100 as Bi->infty. This is a
     * good approximation for this problem since it results in the outside of
     * the can reaching the external temperature incredibly quickly (<1sec). */
    if(Bi<100.00)
        ApplyNaturalBC1D(p, 0, &IsOnRightBoundary, &ConvBC);
    else
        ApplyEssentialBC1D(p, 0, &IsOnRightBoundary, &ExternalTemp);

    return;
}

