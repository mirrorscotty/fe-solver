#include <stdlib.h>

#include "basis.h"
#include "mesh1d.h"
#include "matrix.h"
#include "finite-element1d.h"
#include "solution.h"
#include "mtxsolver.h"

struct fe1d* CreateFE1D(basis *b,
        Mesh1D* mesh,
        matrix* (*makedj)(struct fe1d*, Elem1D*, matrix*),
        matrix* (*makej)(struct fe1d*, Elem1D*, matrix*),
        matrix* (*makef)(struct fe1d*, Elem1D*, matrix*),
        void (*applybcs)(struct fe1d*),
        int maxsteps)
{
    struct fe1d* problem;
    problem = (struct fe1d*) calloc(1, sizeof(struct fe1d));

    problem->b = b;
    problem->mesh = mesh;
    problem->makedj = makedj;
    problem->makej = makej;
    problem->makef = makef;
    problem->applybcs = applybcs;

    problem->tol = 1e-10;

    /* Todo: fix this so it works always, not just for linear interpolation. */
    problem->nrows = (mesh->nelem+1);

    /* Initialize values. These should probably be set later in the program. */
    problem->nvars = 0;

    problem->t = 0;
    problem->maxsteps = maxsteps;
    problem->soln = (solution**) calloc(problem->maxsteps, sizeof(solution*));

    problem->F = NULL;
    problem->R = NULL;

    problem->J = NULL;
    problem->dJ = NULL;

    return problem;
}

void DestroyFE1D(struct fe1d* p)
{
    if(p) {
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
    }
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
        DestroyMatrix(j);
    }

    problem->J = J;

    return J;
}

matrix* AssembledJ1D(struct fe1d *problem, matrix *guess)
{
    Mesh1D *mesh;
    mesh = problem->mesh;

    basis *b;
    b = problem->b;

    matrix *dJ, *dj;
    int i, x, y;
    int c, d;

    double rows = problem->nrows*problem->nvars;

    if(!guess)
        guess = CreateMatrix(rows, 1);

    dJ = CreateMatrix(rows, rows);

    for(i=0; i<mesh->nelem; i++) {
        dj = problem->makedj(problem, mesh->elem[i], guess);

        for(x=0; x<b->n; x++) {
            for(y=0; y<b->n; y++) {
                c = valV(mesh->elem[i]->map, x);
                d = valV(mesh->elem[i]->map, y);

                addval(dJ, val(dj, x, y), c, d);
            }
        }
        DestroyMatrix(dj);
    }

    problem->dJ = dJ;

    return dJ;
}

matrix *AssembleF1D(struct fe1d *problem, matrix *guess)
{
    Mesh1D *mesh = problem->mesh;
    basis *b = problem->b;

    matrix *F, *f;// *lguess;
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

/* Initialize the problem for the transient solver. This function assumes the
 * equations to be solved are linear.
 */
void FE1DTransInit(struct fe1d *problem, matrix* InitSoln)
{
    matrix *f;
    matrix *dsoln;

    /* Store the initial solution in the load vector so that the boundary
     * conditions get applied to it. */
    problem->F = InitSoln;

    AssembleJ1D(problem, NULL);
    AssembledJ1D(problem, NULL);

    /* Apply the boundary conditions */
    problem->applybcs(problem);

    /* Generate the appropriate load vector */
    AssembleF1D(problem, NULL);

    /* Apply the boundary conditions again to make sure the load vector is
     * properly initialized */
    problem->applybcs(problem);

    /* Solve for the time derivative at t=0 */
    f = mtxmul(problem->J, InitSoln);
    mtxneg(f);
    dsoln = SolveMatrixEquation(problem->dJ, f);
    DestroyMatrix(f);

    /* Store the initial solution */
    StoreSolution(problem, InitSoln, dsoln);
}

/*
   matrix *AssembleF1DTrans(struct fe1d *problem, matrix *guess)
   {
   Mesh1D *mesh = problem->mesh;
   basis *b = problem->b;

   matrix *M, *m;
   solution *prevsoln;
   int i, x, y;
   int c, d;

   double rows = problem->nrows*problem->nvars;

   if(!guess)
   guess = CreateMatrix(rows, 1);

   M = CreateMatrix(rows, rows);

   for(i=0; i<mesh->nelem; i++) {
   m = problem->makef(problem, mesh->elem[i], guess);

   for(x=0; x<b->n; x++) {
   for(y=0; y<b->n; y++) {
   c = valV(mesh->elem[i]->map, x);
   d = valV(mesh->elem[i]->map, y);

   addval(M, val(m, x, y), c, d);
   }
   }
   DestroyMatrix(m);
   }

   prevsoln = FetchSolution(problem, problem->t-1);

   problem->F = mtxmul(M, prevsoln->val);

   return M;
   }
   */

/* These functions are copied almost verbatim from the 2d file. */
/* Function to apply Dirchlet boundary conditions. The second argument is a
   function that determines whether or not the boundary condition is applied
   for the current row (node). If the function evaluates as true, then the BC
   is applied. The value at that node is then set to whatever the other
   function returns. */
void ApplyEssentialBC1D(struct fe1d* p,
        int var,
        int (*cond)(struct fe1d*, int),
        double (*BC)(struct fe1d*, int))
{
    int i, j;
    /* Loop through the rows */
    for(i=var; i<p->nrows*p->nvars; i+=p->nvars) {
        /* Check to see if the BC should be applied */
        if(cond(p, i)) {
            //printf("Applying Dirchlet boundary condition at node %d for variable %d\n", i/p->nvars, var);
            /* Zero out the row */
            for(j=0; j<p->nrows*p->nvars; j++) {
                setval(p->J, (i==j)?1:0, i, j); /* Use the Chroniker delta */
                setval(p->dJ, (i==j)?1:0, i, j); /* Use the Chroniker delta */
            }
            /* Set the appropriate value in the load vector */
            setval(p->F, BC(p, i), i, 0);
        }
    }
    return;
}

/* This works the same way, only it applies a Neumann BC */
void ApplyNaturalBC1D(struct fe1d *p,
        int var,
        int (*cond)(struct fe1d*, int),
        double (*BC)(struct fe1d*, int))
{
    int i;
    for(i=var; i<p->nrows*p->nvars; i+=p->nvars) {
        if(cond(p, i))
            //printf("Applying Neumann boundary condition at node %d for variable %d\n", i/p->nvars, var);
            addval(p->F, BC(p, i), i, 0);
    }
}

/* Apply the initial condition for the specifed variable. The supplied function
 * should return the value of the dependant variable and accept the independant
 * variable as its only argument. */
matrix* GenerateInitialCondition(struct fe1d *p,
        int var,
        double (*f)(double))
{
    int i;
    double x; /* The x value at the current mesh node. */
    int n = p->nvars; /* The total number of dependant variables */

    matrix *InitSoln;
    InitSoln = CreateMatrix(p->nrows*p->nvars, 1);

    for(i=var; i<p->nrows*n; i+=n) {
        x = valV(p->mesh->nodes, i/n);
        setval(InitSoln, f(x), i, 0);
    }

    return InitSoln;
    /* This crap here is to make sure that the boundary conditions are applied
     * to the initial state of the system. Bad things happen if this isn't done.
     */
    //p->F = InitSoln;
    //    p->J = CreateMatrix(p->nrows*p->nrows, p->nrows*p->nrows);
    //  p->applybcs(p);
    //DestroyMatrix(p->J);

    //StoreSolution(p, InitSoln);
}

/* Store the solution and advance the current time index. Returns 0 on failure.
 * The first arguement is the problem to store the solutions in.
 * Argument 2: the values to store.
 * Argument 3: The first time derivative. If this is NULL, it is ignored.
 */
int StoreSolution(struct fe1d* p, matrix* values, matrix* dvalues)
{
    if(p->t == p->maxsteps)
        return 0;
    p->soln[p->t] = CreateSolution(p->t, p->dt, values);

    /* Store the time derivative if it is supplied. */
    if(dvalues)
        p->soln[p->t]->dval = dvalues;

    p->t = p->t + 1;
    return 1;
}

solution* FetchSolution(struct fe1d *p, int t)
{
    if(t < 0)
        return NULL;

    if(t < p->t)
        return p->soln[t];
    else
        return NULL;
}

void PrintSolution(struct fe1d *p, int t)
{
    solution *s;
    s = FetchSolution(p, t);
    mtxprnt(s->val);
}

