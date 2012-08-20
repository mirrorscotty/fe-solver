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

#include "material-data/freezing/freezing.h"

/* The following two functions are for heat conduction. */
/* Creates the Jacobian and helps solve for the current time step */
double Residual(struct fe1d *p, matrix *guess, Elem1D *elem, double x, int f1, int f2)
{
    double value = 0;
    basis *b;
    b = p->b;
    
    value  = b->phi[f1](x) * b->phi[f2](x);

    return value;
}

/* Used in calculating the load vector/stuff from the previous time step. */
/* Note: This function fails pretty badly if a material data file isn't loaded.
 */
double R2D2(struct fe1d *p, matrix *guess, Elem1D *elem, double x, int f1, int f2)
{
    double value = 0;
    basis *b;
    b = p->b;

    /* This is just to fetch the value of T */
    double T;
    solution *s;
    s = FetchSolution(p, p->t-1);
    T = EvalSoln1D(p, 0, elem, s, x);

    value  = b->phi[f1](x) * b->phi[f2](x);
    value += p->dt * IMap1D(p, elem, x)
             * b->dphi[f1](x) * b->dphi[f2](x);// * alpha(T);

    return value;
}

/* The next two are for reaction */
double ReactResidual(struct fe1d *p, matrix *guess, Elem1D *elem, double x, int f1, int f2)
{
    return p->b->phi[f1](x) * p->b->phi[f2](x);
}

double R5D2(struct fe1d *p, matrix *guess, Elem1D *elem, double x, int f1, int f2)
{
    double value;
    double R = 8.314, A = 1, Ea = 2000; /* Default values to make gcc happy. */

    /* This is just to fetch the value of T */
    double T;
    solution *s;
    s = FetchSolution(p, p->t-1);
    T = EvalSoln1D(p, 0, elem, s, x);

    value = p->b->phi[f1](x) * p->b->phi[f2](x);
    value += p->dt * A * exp(-Ea / (R * T));

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
            value = quad1d3generic(p, guess, elem, &R2D2, i/v, j/v);
            setval(m, value, i, j);

            //value = quad1d3generic(p, guess, elem, &ReactResidual, i/v, j/v);
            //setval(m, value, i+1, j+1);
        }
    }

    return m;
}

/* Create the load vector... thing */
matrix* CreateElementLoad(struct fe1d *p, Elem1D *elem, matrix *guess) {
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

            //value = quad1d3generic(p, guess, elem, &R5D2, i/v, j/v);
            //setval(m, value, i+1, j+1);
        }
    }

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

double ConvBC(struct fe1d *p, int row)
{
    double h = 1;
    double Tinf = 290;
    double T;

    /* Fetch the value of T from the previous solution. If this is being
     * applied at the initial condition, then simply return 0.*/
    solution *s;
    s = FetchSolution(p, p->t-1);
    if(s) {
        T = val(s->values, row, 0);
        printf("T = %g\n", T);
        return -h*(T-Tinf)/T;
    } else {
        return 0;
    }
}

void ApplyAllBCs(struct fe1d *p)
{
    
    // BC at x=0
    ApplyEssentialBC1D(p, 0, &IsOnLeftBoundary, &Left);
    
    // BC at x=L
    ApplyEssentialBC1D(p, 0, &IsOnRightBoundary, &Right);
    
    return;
}


/* Function that sets the initial condition. */
double InitTemp(double x)
{
    return 273;
}

double InitC(double x)
{
    return 1;
}

int main(int argc, char *argv[])
{
    Mesh1D *mesh;
    basis *b;
    struct fe1d* problem;

    /* Load a data file if one is supplied. */
    if(argc == 2)
        init(argv[1]);

    /* Make a linear 1D basis */
    b = MakeLinBasis(1);

    /* Create a uniform mesh */
    mesh = GenerateUniformMesh1D(b, 0.0, 1.0, 5);
    
    problem = CreateFE1D(b, mesh, &CreateElementMatrix, &CreateElementLoad, &ApplyAllBCs, 100);
    problem->nvars = 1;
    problem->dt = .1;

    ApplyInitialCondition(problem, 0, &InitTemp); /* Initial temperature */
    //ApplyInitialCondition(problem, 1, &InitC); /* Initial concentration */

    while(problem->t<problem->maxsteps)
        StoreSolution(problem, LinSolve1DTrans(problem));
    PrintSolution(problem, problem->maxsteps-1);

    return 0;
}
