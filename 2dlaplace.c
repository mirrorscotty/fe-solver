#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "basis.h"
#include "mtxsolver.h"
#include "integrate.h"
#include "mesh.h"
#include "finite-element.h"

matrix* CreateElementMatrix(struct fe *p, matrix *guess)
{
    basis *b;
    Elem2D *elem;
    
    b = p->b;
    /* Todo: Fix this. */
    elem = &(p->mesh->elem[0]);
    
    int i, j;
    double value = 0;
    double dx = elem->dx;
    double dy = elem->dy;
    matrix *m;

    m = CreateMatrix(b->n, b->n);

    for(i=0; i<b->n; i++) {
        for(j=0; j<b->n; j++) {
            value = dy/dx * quad2d3(b, i, 1, 0) * quad2d3(b, j, 1, 0);
            value += dx/dy * quad2d3(b, i, 0, 1) * quad2d3(b, j, 0, 1);
            setval(m, value, i, j);
        }
    }

    return m;
}

/* Create the load vector */
matrix* CreateElementLoad(struct fe *p, matrix *guess) {
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
void ApplyBC(matrix* J, matrix* F, basis *b, int node, double value)
{
    int rows = mtxlen2(J);
    int i;
    int o = b->overlap;

    for(i=0;i<rows; i++) {
        setval(J, 0, node, i);
        //setval(J, 0, rows-o, i);
    }

    setval(J, 1, node, node);

    setval(F, value, node, 0);
}

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


int main(int argc, char *argv[])
{
    Mesh2D *mesh;
    basis *b;

    matrix *J, *F, *E;
    
    struct fe* problem;

    b = MakeLinBasis(2);

    /* Create a uniform mesh */
    mesh = GenerateUniformMesh2D(0.0, 1.0,
                                 0.0, 1.0,
                                 2, 2);
    
    problem = CreateFE(b, mesh, CreateElementMatrix, CreateElementLoad, ApplyAllBCs);
    
    AssembleJ(problem, NULL);
    AssembleF(problem, NULL);
    
    problem->applybcs(problem);

    E = SolveMatrixEquation(problem->J, problem->F);

    mtxprnt(E);

    printf("p1 = %g\np2 = %g\np3 = %g\np4 = %g\n",
           EvalLin2D(b, 0, 0, 0),
           EvalLin2D(b, 1, 0, 1),
           EvalLin2D(b, 2, 1, 0),
           EvalLin2D(b, 3, 1, 1));

    //J = AssembleJ(&CreateElementMatrix, b, mesh, Pe);
    //J = AssembleJ(&testelem, b, mesh, Pe);

    //F = AssembleF(&CreateElementLoad, b, mesh, Pe);
    //F = AssembleF(&testload, b, mesh, Pe);

    //ApplyBoundaryConditions(J, F, b, 0, 1);
    //ApplyBoundaryConditions(J, F, b, 2, 1);
    
    //E = SolveMatrixEquation(J, F);

    /* Clean up the allocated memory */
    //DestroyMatrix(J);
    //DestroyMatrix(F);
    //DestroyMatrix(E);

    return 0;
}
