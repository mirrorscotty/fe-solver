#include <stdlib.h>

#include "basis.h"
#include "mesh.h"
#include "matrix.h"
#include "finite-element.h"

struct fe* CreateFE(basis* b,
                    Mesh2D* mesh,
                    matrix* (*makej)(struct fe*, Elem2D*, matrix*),
                    matrix* (*makef)(struct fe*, Elem2D*, matrix*),
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
 * argument is a structure that describes the FE problem, and the second
 * is the initial guess.
 */
matrix* AssembleJ(struct fe *problem, matrix *guess)
{
    Mesh2D *mesh;
    mesh = problem->mesh;

    basis *b;
    b = problem->b;

    matrix *J, *j;
    int i, x, y;
    int c, d;
    //int z = 0;

    /* Determine the number of rows in the coefficient matrix. This expression
     * needs to change to be able to accomodate biquad elements. */
    double rows = (mesh->nelemx+1)*(mesh->nelemy+1)*problem->nvars;
   
    /* If an initial guess is not supplied, then the initial guess is all
     * zeroes */
    if(!guess)
        guess = CreateMatrix(rows, 1);

    /* Create the blank coefficient matrix */
    J = CreateMatrix(rows, rows);

    /* Iterate through the list of elements in the mesh. Create the element
     * matrix for each and then add it to the global matrix as appropriate */
    for(i=0; i<mesh->nelemx*mesh->nelemy; i++) {
        /* Generate the element matrix for the specified element width */
    //MeshPrint(mesh);
        j = problem->makej(problem, mesh->elem[i], guess);

        //j = CreateOnesMatrix(8, 8); // for testing purposes

        /* Add the values of the element matrix to the global matrix */
        for(x=0; x<b->n; x++) {
            for(y=0; y<b->n; y++) {
                c = valV(mesh->elem[i]->map, x);
                d = valV(mesh->elem[i]->map, y);

                addval(J, val(j, x, y), c, d);
                //addval(J, val(j, x*2, y*2), c*2, d*2);
                //addval(J, val(j, x*2+1, y*2), c*2+1, d*2);
                //addval(J, val(j, x*2, y*2+1), c*2, d*2+1);
                //addval(J, val(j, x*2+1, y*2+1), c*2+1, d*2+1);
            }
        }

        /* Clean up */
        DestroyMatrix(j);
    }
    //J = problem->makej(problem, guess);
    problem->J = J;

    return J;
}

/* Assemble the load vector */
matrix* AssembleF(struct fe *problem, matrix *guess)
{
    Mesh2D *mesh = problem->mesh;
    basis *b = problem->b;
    
    matrix *F, *f, *lguess;
    int i, j, rows;
    int r = b->n - b->overlap;

    /* Determine the number of rows required for the global matrix */
    rows = (mesh->nelemx+1)*(mesh->nelemy+1);

    F = CreateMatrix(rows, 1);

    for(i=0; i<rows-r; i=i+r) {
        lguess = GetLocalGuess(problem, guess, i);
        f = problem->makef(problem, mesh->elem[i], lguess);
        DestroyMatrix(lguess);

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

matrix* GetLocalGuess(struct fe *p, matrix *guess, int elem)
{
    matrix* lguess;
    int x, y, i, j;
    lguess = CreateMatrix(p->nvars*p->b->n, 1);
    
    int ny = p->mesh->nelemy;
    
    for(i=0; i<mtxlen2(lguess); i+=p->nvars) {
        x = elem / (p->mesh->nelemy+1);
        y = elem - x*(p->mesh->nelemy+1);
        
        for(j=0; j<p->nvars; j++) {
            setval(lguess, val(guess, p->nvars*(x)+j, 0), 0, 0);
            setval(lguess, val(guess, p->nvars*(x+1)+j, 0), p->nvars+j, 0);
            setval(lguess, val(guess, p->nvars*(x+1+ny)+j, 0), 2*p->nvars+j, 0);
            setval(lguess, val(guess, p->nvars*(x+1+ny+1)+j, 0), 3*p->nvars+j, 0);
        }
    }

    return lguess;
}

/* Function to apply Dirchlet boundary conditions. The second argument is a
   function that determines whether or not the boundary condition is applied
   for the current row (node). If the function evaluates as true, then the BC
   is applied. The value at that node is then set to whatever the other
   function returns. */
void ApplyEssentialBC(struct fe* p,
                      int (*cond)(struct fe*, int),
                      double (*BC)(struct fe*, int))
{
    int i, j;
    /* Loop through the rows */
    for(i=0; i<mtxlen2(p->J); i++) {
        /* Check to see if the BC should be applied */
        if(cond(p, i)) {
            /* Zero out the row */
            for(j=0; j<mtxlen2(p->J); j++) 
                setval(p->J, (i==j)?1:0, i, j); /* Use the Chroniker delta */
            /* Set the appropriate value in the load vector */
            setval(p->F, BC(p, i), i, 0);
        }
    }
    return;
}
