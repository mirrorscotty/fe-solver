#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "basis.h"
#include "integrate.h"
#include "mesh2d.h"
#include "finite-element.h"

matrix* CreateElementMatrix(struct fe *p, Elem2D* e, matrix *guess)
{
    basis *b;
    
    b = p->b;
    
    int i, j;
    double value = 0;
    matrix *m;

    m = CreateMatrix(b->n, b->n);

    for(i=0; i<b->n; i++) {
        for(j=0; j<b->n; j++) {
            value = quad2d32d3(p, e, i, j, 1, 0);
            value += quad2d32d3(p, e, i, j, 0, 1);
            setval(m, value, i, j);
        }
    }

    return m;
}

/* Create the load vector */
matrix* CreateElementLoad(struct fe *p, Elem2D *e, matrix *guess) {
    int i;
    matrix *f;
    
    basis *b;
    b = p->b;

    e = NULL; /* Just to make the compiler stop complaining about unused vars */

    f = CreateMatrix(b->n, 1);

    for(i=0; i<b->n; i++) {
        setval(f, 0, i, 0);
    }

    return f;
}

int IsOnInnerRadius(struct fe* p, int node)
{
    double r, theta, x, y;
    double f, a, b, e;
    double tol = 1e-5;
    vector *v;

    e = 0.8; /* Eccentricity */

    b = 1;
    a = sqrt((1-e*e)*b*b);

    v = p->mesh->nodes[node];
    x = valV(v, 0);
    y = valV(v, 1);

    theta = atan(y/x);
    r = sqrt(x*x+y*y);

    f = 1/sqrt(pow(cos(theta), 2)/(a*a) + pow(sin(theta), 2)/(b*b));

    if(fabs(f-r) < tol)
        return 1;
    else
        return 0;
}

int IsOnOuterRadius(struct fe* p, int node)
{
    double r, x, y;
    double tol = 1e-5;
    vector *v;

    double f = 5; /* r* */

    v = p->mesh->nodes[node];

    x = valV(v, 0);
    y = valV(v, 1);

    r = sqrt(x*x + y*y);

    if(fabs(f-r) < tol) {
        return 1;
    }
    else
        return 0;
}

double InnerRadiusBC(struct fe *p, int row)
{
    return 1;
}

double OuterRadiusBC(struct fe *p, int row)
{
    double x = valV(p->mesh->nodes[row], 0);
    double y = valV(p->mesh->nodes[row], 1);

    return 1/sqrt(x*x+y*y);
}

void ApplyAllBCs(struct fe *p)
{
    /* BC at r = the surface */
    ApplyEssentialBC(p, &IsOnInnerRadius, &InnerRadiusBC);

    /* BC at r = r* */
    ApplyEssentialBC(p, &IsOnOuterRadius, &OuterRadiusBC);
}

void PrintResults(struct fe *p, matrix *soln)
{
    matrix* tmp;
    int i;
    tmp = CreateMatrix(nRows(soln), 3);
    for(i=0; i<nRows(soln); i++) {
        setval(tmp, valV(p->mesh->nodes[i], 0), i, 0);
        setval(tmp, valV(p->mesh->nodes[i], 1), i, 1);
        setval(tmp, val(soln, i, 0), i, 2);
    }
    mtxprnt(tmp);
    DestroyMatrix(tmp);
    return;
}

int main(int argc, char *argv[])
{
    Mesh2D *mesh;
    basis *b;

    matrix *E;
    
    struct fe* problem;

    b = MakeQuadBasis(2);

    /* Create a uniform mesh */
    mesh = MakeSpheroidMesh(b, 0.8, 5, 10, 10);
    
    problem = CreateFE(b, mesh, CreateElementMatrix, CreateElementLoad, ApplyAllBCs);
    problem->nvars = 1;
    
    E = LinSolve(problem);

//    mtxprnt(E);
    PrintResults(problem, E);

    DestroyFE(problem);
    //DestroyMatrix(E);

    return 0;
}
