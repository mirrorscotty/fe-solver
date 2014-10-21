#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "basis.h"
#include "integrate.h"
#include "mesh2d.h"
#include "finite-element.h"

matrix* CreateElementMatrix(struct fe *p, Elem2D *elem, matrix *guess)
{
    basis *b;
    
    b = p->b;
    
    int i, j;
    double value = 0;
    double dx = 1;//elem->dx;
    double dy = 1;//elem->dy;
    matrix *m;

    m = CreateMatrix(b->n, b->n);

    for(i=0; i<b->n; i++) {
        for(j=0; j<b->n; j++) {
            value = dy/dx * quad2d32d3(p, elem, i, j, 1, 0);
            value += dx/dy * quad2d32d3(p, elem, i, j, 0, 1);
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

/* Simple function to alter J and F so that Dirichlet boundary conditions are
 * imposed on both ends of the domain. In this case, "leftbc" is imposed at
 * c = 0, and rightbc is imposed at c = 1.
 */
/*
void ApplyBC(matrix* J, matrix* F, basis *b, int node, double value)
{
    int rows = nRows(J);
    int i;
    int o = b->overlap;

    for(i=0;i<rows; i++) {
        setval(J, 0, node, i);
        //setval(J, 0, rows-o, i);
    }

    setval(J, 1, node, node);

    setval(F, value, node, 0);
}
*/
/*
void ApplyAllBCs(struct fe *p)
{
    matrix *J, *F;
    basis *b;
    Mesh2D *mesh;
    
    b = p->b;
    mesh = p->mesh;
    F = p->F;
    J = p->J;
    
    int i;
    vector *v;

    // BC at y=0
    for(i=0; i<mesh->nelemx+1; i++) {
        v = GetNodeCoordinates(mesh, i);
        ApplyBC(J, F, b, i, valV(v,1)*(2-valV(v,1)) );
        DestroyVector(v);
    }

    // BC at x=0
    for(i=0; i< (mesh->nelemx+1)*(mesh->nelemy+1); i+= (mesh->nelemy+1)) {
        v = GetNodeCoordinates(mesh, i);
        ApplyBC(J, F, b, i, valV(v,0)*(valV(v,0)-2) );
        DestroyVector(v);
    }
}
*/
int IsOnXAxis(struct fe* p, int node)
{
    double tol = 1e-10;
    if(fabs(valV(p->mesh->nodes[node], 1)) < tol)
        return 1;
    else
        return 0;
}

int IsOnYAxis(struct fe* p, int node)
{
    double tol = 1e-10;
    if(fabs(valV(p->mesh->nodes[node], 0)) < tol)
        return 1;
    else
        return 0;
}

double YAxisBC(struct fe *p, int node)
{
    double y = valV(p->mesh->nodes[node], 1);
    return y*(2-y);
}

double XAxisBC(struct fe *p, int node)
{
    double x = valV(p->mesh->nodes[node], 0);
    return x*(x-2);
}

void ApplyAllBCs2(struct fe *p)
{
    ApplyEssentialBC(p, &IsOnXAxis, &XAxisBC);
    ApplyEssentialBC(p, &IsOnYAxis, &YAxisBC);
}

int main(int argc, char *argv[])
{
    Mesh2D *mesh;
    basis *b;

    matrix *E;
    
    struct fe* problem;

    b = MakeLinBasis(2);

    /* Create a uniform mesh */
    mesh = GenerateUniformMesh2D(0.0, 1.0,
                                 0.0, 1.0,
                                 1, 1);
    
    MeshPrint(mesh);
    problem = CreateFE(b, mesh, CreateElementMatrix, CreateElementLoad, ApplyAllBCs2);
    problem->nvars = 1;
    
    //E = SolveMatrixEquation(problem->J, problem->F);
    E = LinSolve(problem);

    mtxprnt(E);

    printf("p1 = %g\np2 = %g\np3 = %g\np4 = %g\n",
           EvalLin2D(b, 0, 0, 0),
           EvalLin2D(b, 1, 0, 1),
           EvalLin2D(b, 2, 1, 0),
           EvalLin2D(b, 3, 1, 1));

    return 0;
}
