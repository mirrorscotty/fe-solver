#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "basis.h"
#include "mtxsolver.h"
#include "integrate.h"
#include "mesh1d.h"
#include "isoparam.h"
#include "finite-element1d.h"
#include "auxsoln.h"

#include "material-data/freezing/freezing.h"

#include "heat-gui.h"

extern double To;

/* The following two functions are for heat conduction. */
/* Creates the Jacobian and helps solve for the current time step */
double Residual(struct fe1d *p, matrix *guess, Elem1D *elem, double x, int f1, int f2)
{
    double value = 0;
    basis *b;
    b = p->b;

     /* This is just to fetch the value of T */
    //double T;
    //solution *s;
    //s = FetchSolution(p, p->t-1);
    //T = EvalSoln1D(p, 0, elem, s, x);

   
    value  = b->dphi[f1](x) * b->dphi[f2](x);
//    value *= 1e-5;
    value *= IMap1D(p, elem, x);
    value *= IMapCyl1D(p, elem, x);

//    if(valV(elem->map, f1) == p->mesh->nnodes-1) {
//        value += ConvBC(p, valV(elem->map, f1));
//    }

    return value;
}

/* Used in calculating the load vector/stuff from the previous time step. */
/* Note: This function fails pretty badly if a material data file isn't loaded.
 */
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

int IsOnRightBoundary(struct fe1d *p, int row)
{
    double width = p->mesh->x2 - p->mesh->x1;
    double x = valV(p->mesh->nodes, row/p->nvars);
    
    if(fabs(x - width) < 1e-5)
        return 1;
    else
        return 0;
}

int IsOnLeftBoundary(struct fe1d *p, int row)
{
    double x = valV(p->mesh->nodes, row/p->nvars);
  
    if(fabs(x) < 1e-5)
        return 1;
    else
        return 0;
}

double Left(struct fe1d *p, int row)
{
    return 273.0;
}

double Right(struct fe1d *p, int row)
{
    return 290.0;
}

double Zero(struct fe1d *p, int row)
{
    return 0.0;
}

/* The way this function is implemented is probably not mathematically accurate.
 * Strictly speaking, for the implicit solver, the temperature fetched should be
 * the temperature at the next time step, not the one that we already have the
 * solution for. The way it is now should be good enough (tm).*/
double ConvBC(struct fe1d *p, int row)
{
    double h = -1;
    double Tinf = 274;
    double T;

    /* Fetch the value of T from the previous solution. If this is being
     * applied at the initial condition, then simply return 0.*/
    solution *s;
    s = FetchSolution(p, p->t-1);
    if(s) {
        T = val(s->val, row, 0);
        printf("T = %g\n", T);
        if(T==0)
            return 0;
        return h*(T-Tinf);
    } else {
        return 0;
    }
}

void ApplyAllBCs(struct fe1d *p)
{
    
    // BC at x=0
    //ApplyEssentialBC1D(p, 0, &IsOnLeftBoundary, &Left);
    
    // BC at x=L
    ApplyNaturalBC1D(p, 0, &IsOnRightBoundary, &ConvBC);

    return;
}

/* ODEs */
double react1(double cprev, double T, double dt)
{
    double A, Ea, R;
    A = 1;
    Ea = 1;
    R = 1;

    T = fabs(T);
    
    return cprev - dt*A*exp(-Ea/(R*T));
}

double react2(double cprev, double T, double dt)
{
    double A, Ea, R;
    A = 2;
    Ea = 2;
    R = 2;

    T = fabs(T);
    
    return cprev - dt*A*exp(-Ea/(R*T));
}

