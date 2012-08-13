#include <stdlib.h>

#include "basis.h"
#include "mesh1d.h"
#include "matrix.h"
#include "finite-element1d.h"

struct fe1d* CreateFE1D(basis *b,
                        Mesh1D* mesh,
                        matrix* (*makej)(struct fe1d*, Elem1D*, matrix*),
                        matrix* (*makef)(struct fe1d*, Elem1D*, matrix*),
                        void (*applybcs)(struct fe1d*))
    {
    struct fe1d* problem;
    problem = (struct fe1d*) calloc(1, sizeof(struct fe1d));

    problem->b = b;
    problem->mesh = mesh;
    problem->makej = makej;
    problem->makef = makef;
    problem->applybcs = applybcs;

    problem->tol = 1e-10;

    /* Todo: fix this so it works always, not just for linear interpolation. */
    problem->nrows = (mesh->nelem+1);

    /* Initialize values. These should probably be set later in the program. */
    problem->nvars = 0;

    problem->F = NULL;
    problem->R = NULL;
    problem->J = NULL;

    return problem;
}

void DestroyFE1D(struct fe1d* p)
{
    if(p->b)
        DestroyBasis(p->b);
    if(p->F)
        DestroyMatrix(p->F);
    if(p->R)
        DestroyMatrix(p->R);
    if(p->J)
        DestroyMatrix(p->J);
    if(p->mesh)
        DestroyMesh1D(p->mesh);
    return;
}

matrix* AssembleJ1D(struct fe1d *problem, matrix *guess)
{
    Mesh1D *mesh;
    mesh = problem->mesh;

    basis *b;
    b = problem->b;

    matrix *J, *j;
    int i, x, y;
    int c, d;

    double rows = problem->nrows*problem->nvars;

    if(!guess)
        guess = CreateMatrix(rows, 1);

    J = CreateMatrix(rows, rows);

    for(i=0; i<mesh->nelem; i++) {
        j = problem->makej(problem, mesh->elem[i], guess);

        for(x=0; x<b->n; x++) {
            for(y=0; y<b->n; y++) {
                c = valV(mesh->elem[i]->map, x);
                d = valV(mesh->elem[i]->map, y);

                addval(J, val(j, x, y), c, d);
            }
        }
    }
    DestroyMatrix(j);

    problem->J = J;

    return J;
}

matrix *AssembleF1D(struct fe1d *problem, matrix *guess)
{
    Mesh1D *mesh = problem->mesh;
    basis *b = problem->b;

    matrix *F, *f, *lguess;
    int i, j, rows;
    int r = b->n - b->overlap;

    rows = problem->nrows;

    F = CreateMatrix(rows, 1);

    for(i=0; i<rows-r; i=i+r) {

        f = problem->makef(problem, mesh->elem[i], guess);
        
        for(j=0; j<b->n; j++) {
            addval(F, val(f, j, 0), i+j, 0);
        }

        DestroyMatrix(f);
    }
    problem->F = F;

    return F;
}

