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
    double h = 1e-5;
    double DkDx = 0, Ti, ki;
    int i;
    basis *b;
    b = p->b;

    /* This is just to fetch the value of T */
    double T, Tu;

    solution *s;
    s = CreateSolution(p->t, p->dt, guess);
    T = EvalSoln1D(p, 0, elem, s, x);

    /* Get the real value of T, not the dimensionless temperature. */
    Tu = uscaleTemp(p->charvals, T);

    term1 = k(comp_global, Tu)/p->charvals.k;// * b->phi[f1](x);
    term1 *= b->dphi[f1](x) * b->dphi[f2](x);
    term1 *= IMap1D(p, elem, x);

    /* Calculate gradient of thermal conductivity. Hopefully this is the
     * mathematically correct way to handle this. */
    for(i=0; i<b->n; i++) {
        Ti = EvalSoln1D(p, 0, elem, s, valV(elem->points, i));
        ki = k(comp_global, uscaleTemp(p->charvals, Ti));
        DkDx += ki * b->dphi[f1](x);
    }
    /* Now that we've calculated T, we no longer need this */
    free(s);

//    term2 = (k(comp_global, Tu+h)-k(comp_global, Tu-h))/(2*h);
//    term2 *= 1/p->charvals.k;
//    term2 *= pow(b->dphi[f1](x), 2) * b->phi[f2](x) * T;
    term2 = DkDx/p->charvals.k * b->dphi[f1](x) * b->phi[f2](x);
    term2 /= IMap1D(p, elem, x);
    /* Hopefully the above term is correct. It appears to yield accurate
     * results, at least. */
    
    return term2-term1;
    //return -1*term1;
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
            //* (p->charvals.alpha/p->charvals.k);
            / p->charvals.alpha;
    value *= 1/IMap1D(p, elem, x);

    return -1*value;
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
    /*
    double width = p->mesh->x2 - p->mesh->x1;
    double x = valV(p->mesh->nodes, row/p->nvars);
    
    if(fabs(x - width) < 1e-7)
        return 1;
    else
        return 0;
    */
    if(row == len(p->mesh->nodes)-1)
        return 1;
    else
        return 0;
}

/* Same as above, only for the left boundary */
int IsOnLeftBoundary(struct fe1d *p, int row)
{
    /*
    double x = valV(p->mesh->nodes, row/p->nvars);
  
    if(fabs(x) < 1e-7)
        return 1;
    else
        return 0;
    */
    if(row == 0)
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
        //if(T==0)
        //    return 0;
        //else
            return Bi*(T-Tinf);
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

/* Calculate the deformation gradient of the sample at (x,t) based on density
 */
double DeformationGrad(struct fe1d *p, double x, double t)
{
    double T0, Tn, rho0, rhon;
    solution *s0, *sn;
    
    s0 = FetchSolution(p, 0);
    sn = FetchSolution(p, t);
    Tn = uscaleTemp(p->charvals, EvalSoln1DG(p, 0, sn, x, 0));
    T0 = uscaleTemp(p->charvals, EvalSoln1DG(p, 0, s0, x, 0));
    rho0 = rho(comp_global, T0);
    rhon = rho(comp_global, Tn);

    return rhon/rho0;
}

