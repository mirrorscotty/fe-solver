/**
 * @file heat-gui.c
 * Set of functions to solve the heat equation with reactions occuring.
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
#include "auxsoln.h"

#include "material-data/choi-okos/choi-okos.h"

#include "heat-gui.h"

extern double EaA, EaB, AA, AB,
       Text_hot, Text_cold, t_heat;
extern choi_okos *comp_global;
//double h = 5; // Not actually used.

/* Take the values stored in global variables for heating time and temperatures
 * and construct a function for the external temperature as a function of time.
 */
double T_ext(double time)
{
    if(time > t_heat)
        return Text_cold;
    else
        return Text_hot;
}

/* The following two functions are for heat conduction. */
/* Creates the Jacobian and helps solve for the current time step */
double Residual(struct fe1d *p, matrix *guess, Elem1D *elem, double x, int f1, int f2)
{
    double value = 0;
    basis *b;
    b = p->b;

    /* This is just to fetch the value of T */
    double T;
    solution *s;
    //s = FetchSolution(p, p->t-1);
    s = CreateSolution(p->t, p->dt, guess);
    //if(!s)
    if(p->t==0)
        T = p->charvals.Tc;
    else
        T = EvalSoln1D(p, 0, elem, s, x);
    /* Now that we've calculated T, we no longer need this */
    free(s);

    value  = b->dphi[f1](x) * b->dphi[f2](x);
    value *= IMap1D(p, elem, x);
    value *= (alpha(comp_global, uscaleTemp(p->charvals, T))/p->charvals.alpha)
                * b->phi[f1](x) * IMap1D(p, elem, x);
    value *= (alpha(comp_global, uscaleTemp(p->charvals, T))/p->charvals.alpha);
    value *= IMapCyl1D(p, elem, x);

    //printf("alpha_c = %g, alpha = %g, ", p->charvals.alpha, alpha(T));
    //printf("alpha_c/alpha = %g\n", p->charvals.alpha/alpha(T));

    return value;
}

/* Calculate the coefficient matrix for the time derivative unknowns */
double ResDt(struct fe1d *p, matrix *guess, Elem1D *elem, double x, int f1, int f2)
{
    double value;
    basis *b;
    b = p->b;

    value = b->phi[f1](x) * b->phi[f2](x) * 1/IMap1D(p, elem, x);
    value *= IMapCyl1D(p, elem, x);
//    value *= -1;

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
    return scaleTemp(p->charvals, T_ext(uscaleTime(p->charvals, p->t*p->dt)));
}


/* The way this function is implemented is probably not mathematically accurate.
 * Strictly speaking, for the implicit solver, the temperature fetched should be
 * the temperature at the next time step, not the one that we already have the
 * solution for. The way it is now should be good enough (tm).*/
double ConvBC(struct fe1d *p, int row)
{
    double Tinf = scaleTemp(p->charvals, T_ext(uscaleTime(p->charvals, p->dt*p->t)));
    double T;
    double Bi = BiotNumber(p->charvals);

    /* Fetch the value of T from the previous solution. If this is being
     * applied at the initial condition, then simply return 0.*/
    //solution *s;
    //s = FetchSolution(p, p->t-1);
    //if(s) {
    if(p->guess) {
        //T = val(s->val, row, 0);
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

/* ODEs */
double react1(double cprev, double T, double dt)
{
    //double AA, EaA, R;
    double R = 8.314;
    //AA = 1;
    //EaA = 1;

    T = fabs(T);
    
    return cprev*(1 - dt*AA*exp(-EaA/(R*T)));
}

double react2(double cprev, double T, double dt)
{
    //double AB, EaB, R;
    double R = 8.314;
    //AB = 2;
    //EaB = 2;

    T = fabs(T);
    
    return cprev*(1 - dt*AB*exp(-EaB/(R*T)));
}

