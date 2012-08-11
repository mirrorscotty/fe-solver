#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "basis.h"
#include "mtxsolver.h"
#include "integrate.h"
#include "mesh.h"
#include "isoparam.h"
#include "finite-element.h"

double ElemJdRxdx(struct fe *p, matrix *guess, Elem2D *elem, double x, double y, int f1, int f2)
{
    double value;
    
    return value;
}

matrix* CreateElementMatrix(struct fe *p, Elem2D *elem, matrix *guess)
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
            /* dRx/dx */
            value = quad2d3generic(p, guess, elem, &ElemJdRxdx, i/2, j/2);
            setval(m, value, i, j);
        }
    }

    return m;
}

/* Create the load vector */
matrix* CreateElementLoad(struct fe *p, Elem2D *elem, matrix *guess) {
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
int IsOnRightBoundary(struct fe *p, int row)
{
    double width = 1;
    double x = valV(p->mesh->nodes[row/p->nvars], 0);
    
    if(fabs(x - width) < 1e-5)
        return 1;
    else
        return 0;
}

int IsOnLeftBoundary(struct fe *p, int row)
{
    double height = p->mesh->y2 - p->mesh->y1;
    double y = valV(p->mesh->nodes[row/p->nvars], 0);
    
    if(fabs(y - height) < 1e-5)
        return 1;
    else
        return 0;
}

double Zero(struct fe *p, int row)
{
    return 0.0;
}

void ApplyAllBCs(struct fe *p)
{
    
    // BC at x=0
    ApplyNaturalBC(p, 0, &IsOnLeftBoundary, &Zero);
    
    // BC at x=L
    ApplyNaturalBC(p, 0, &IsOnRightBoundary, &Zero);
    
    return;
}

int main(int argc, char *argv[])
{
    Mesh2D *mesh;
    basis *b;
    matrix *E;
    struct fe* problem;

    /* Make a linear 2D basis */
    b = MakeLinBasis(2);

    /* Create a uniform mesh */
    mesh = GenerateUniformMesh1D(0, 1, 5);
    meshprnt(mesh);
    
    //problem = CreateFE(b, mesh, &CreateElementMatrix, &CreateElementLoad, &ApplyAllBCs);
    //problem->nvars = 1;
    
    /* Not used */
    //problem->P = .1;
    //problem->a = .01;
    
    //E = NLinSolve(problem, NULL);

    //mtxprnt(E);
    //puts("");
    
    return 0;
}
