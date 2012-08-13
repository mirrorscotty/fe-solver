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

double Residual(struct fe1d *p, matrix *guess, Elem1D *elem, double x, int f1, int f2)
{
    double value = 0;
    basis *b;
    b = p->b;
    
    value  = IMap1D(p, elem, x)*b->dphi[f1](x);
    value *= b->dphi[f2](x);

    return value;
}

matrix* CreateElementMatrix(struct fe1d *p, Elem1D *elem, matrix *guess)
{
    basis *b;
    b = p->b;
    
    int nvars = p->nvars;
    
    int i, j;
    double value = 0;
    matrix *m;
    
    m = CreateMatrix(b->n*nvars, b->n*nvars);
    
    for(i=0; i<b->n*nvars; i+=nvars) {
        for(j=0; j<b->n*nvars; j+=nvars) {
            value = quad1d3generic(p, guess, elem, &Residual, i, j);
            setval(m, value, i, j);
        }
    }

    return m;
}

/* Create the load vector */
matrix* CreateElementLoad(struct fe1d *p, Elem1D *elem, matrix *guess) {
    int i;
    matrix *f;
    
    basis *b;
    b = p->b;

    f = CreateMatrix(b->n, 1);

    for(i=0; i<b->n; i++) {
        setval(f, 0, i, 0);
    }

    return f;
}

/* Todo: double-check these functions */
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

double Zero(struct fe1d *p, int row)
{
    return 0.0;
}

double One(struct fe1d *p, int row)
{
    return 1.0;
}

void ApplyAllBCs(struct fe1d *p)
{
    
    // BC at x=0
    ApplyEssentialBC1D(p, 0, &IsOnLeftBoundary, &Zero);
    
    // BC at x=L
    ApplyEssentialBC1D(p, 0, &IsOnRightBoundary, &One);
    
    return;
}

int main(int argc, char *argv[])
{
    Mesh1D *mesh;
    basis *b;
    matrix *E;
    struct fe1d* problem;

    /* Make a linear 1D basis */
    b = MakeLinBasis(1);

    /* Create a uniform mesh */
    mesh = GenerateUniformMesh1D(b, 0.0, 1.0, 5);
    
    problem = CreateFE1D(b, mesh, &CreateElementMatrix, &CreateElementLoad, &ApplyAllBCs);
    problem->nvars = 1;
    
    E = LinSolve1D(problem);

    mtxprnt(E);
    
    return 0;
}
