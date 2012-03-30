#include <stdlib.h>

#include "basis.h"
#include "mesh.h"
#include "matrix.h"
#include "finite-element.h"

struct fe* CreateFE(basis* b,
                    Mesh2D* mesh,
                    matrix* (*makej)(struct fe*, matrix*),
                    matrix* (makef)(struct fe*, matrix*),
                    void (*applybcs)(struct fe*))
{
    struct fe* problem;
    problem = (struct fe*) calloc(1, sizeof(struct fe));
    
    problem->b = b;
    problem->mesh = mesh;
    problem->makej = makej;
    problem->makef = makef;
    problem->applybcs = applybcs;
    
    problem->tol = 1e-10;
    
    return problem;
}

void DestroyFE(struct fe* p)
{
    if(p->b)
        DestroyBasis(p->b);
    if(p->F)
        DestroyMatrix(p->F);
    if(p->J)
        DestroyMatrix(p->J);
    if(p->mesh)
        DestroyMesh2D(p->mesh);
    return;
}

/* Assemble the global coefficient matrix for a non-uniform mesh. The first
 * argument is a function pointer to the function which creates the element
 * matrix, and the second is the mesh to be used.
 */
matrix* AssembleJ(struct fe *problem, matrix *guess)
{
    Mesh2D *mesh = problem->mesh;
    basis *b = problem->b;
    matrix *J, *j;
    int nx = mesh->nelemx;
    int ny = mesh->nelemy-1;
    int i, x, y;
    int c, d;

    double rows = (mesh->nelemx+1)*(mesh->nelemy+1);

    J = CreateMatrix(rows, rows);

    for(i=0; i<mesh->nelemx*mesh->nelemy; i++) {
        /* Generate the element matrix for the specified element width */
        j = problem->makej(problem, guess);

        /* Add the values of the element matrix to the global matrix */
        for(x=0; x<b->n; x++) {
            for(y=0; y<b->n; y++) {
                c = i+x+((x>1)?ny:0);
                d = i+y+((y>1)?ny:0);

                c += i/nx;
                d += i/nx;

                addval(J, val(j, x, y), c, d);
            }
        }

        /* Clean up */
        DestroyMatrix(j);
    }
    
    problem->J = J;

    return J;
}

matrix* AssembleF(struct fe *problem, matrix *guess)
{
    Mesh2D *mesh = problem->mesh;
    basis *b = problem->b;
    
    matrix *F, *f;
    int n, i, j, rows;
    int r = b->n - b->overlap;

    /* Determine the number of rows required for the global matrix */
    rows = (mesh->nelemx+1)*(mesh->nelemy+1);

    F = CreateMatrix(rows, 1);

    for(i=0; i<rows-r; i=i+r) {
        f = problem->makef(problem, guess);

        for(j=0; j<b->n; j++) {
            addval(F, val(f, j, 0), i+j, 0);
        }

        DestroyMatrix(f);
    }
    
    problem->F = F;

    return F;
}

matrix* CalcResidual(struct fe *problem, matrix* guess)
{
    matrix *tmp;
    tmp = mtxmul(problem->J, guess);
    problem->R = mtxadd(mtxneg(tmp), problem->F);
    
    
    return problem->R;
}