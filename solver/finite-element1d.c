#include <stdlib.h>
#include <math.h>

#include "basis.h"
#include "mesh1d.h"
#include "matrix.h"
#include "finite-element1d.h"
#include "solution.h"

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

    problem->guess = NULL;

    problem->J = NULL;
    problem->dJ = NULL;

    return problem;
}

void DestroyFE1D(struct fe1d* p)
{
    int i, j;

    if(p) {
        if(p->b)
            DestroyBasis(p->b);
        if(p->F)
            DestroyMatrix(p->F);
        if(p->R)
            DestroyMatrix(p->R);
        if(p->J)
            DestroyMatrix(p->J);
        if(p->dJ)
            DestroyMatrix(p->dJ);
        if(p->mesh)
            DestroyMesh1D(p->mesh);

        if(p->soln) {
            for(i=0; i<p->maxsteps; i++)
                DestroySolution(p->soln[i]);
            free(p->soln);
        }

        if(p->auxsolns) {
            for(i=0; i<p->extravars; i++) {
                for(j=0; j<p->maxsteps; j++)
                    DestroySolution(p->auxsolns[i][j]);
                free(p->auxsolns[i]);
            }
            free(p->auxsolns);
        }

        free(p);
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

    //DestroyMatrix(guess);

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

    //DestroyMatrix(guess);

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

matrix* CalcResidual1D(struct fe1d *problem, matrix* guess)
{
    matrix *tmp;
    tmp = mtxmul(problem->J, guess);
    problem->R = mtxadd(mtxneg(tmp), problem->F);

    return problem->R;
}

/* Initialize the problem for the transient solver. This function assumes the
 * equations to be solved are linear.
 */
void FE1DTransInit(struct fe1d *problem, matrix* InitSoln)
{
    matrix *f;
    matrix *tmp;
    matrix *dsoln;

    /* Store the initial solution in the load vector so that the boundary
     * conditions get applied to it. */
    //problem->F = InitSoln;

    AssembleJ1D(problem, InitSoln);
    AssembledJ1D(problem, InitSoln);

    /* Apply the boundary conditions */
    //problem->applybcs(problem);

    /* Generate the appropriate load vector */
    AssembleF1D(problem, InitSoln);

    /* Apply the boundary conditions again to make sure the load vector is
     * properly initialized */
    problem->applybcs(problem);

    /* Solve for the time derivative at t=0 */
    f = mtxmul(problem->J, InitSoln);
    mtxneg(f);
    tmp = mtxadd(f, problem->F);
    dsoln = SolveMatrixEquation(problem->dJ, tmp);
    DestroyMatrix(f);
    DestroyMatrix(tmp);

    /* Store the initial solution */
    StoreSolution(problem, InitSoln, dsoln);
}

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
}

/**
 * Create a matrix of values to use as the initial condition when solving a
 * transient problem.
 * TODO: Modify this function to allow setting initial conditions for more than
 *      one variable.
 * @param p Finite element structure
 * @param var Number of the variable to set
 * @param value Value to set the variable to at t=0
 */
matrix* GenerateInitCondConst(struct fe1d *p, int var, double value)
{
    int i;
    //double x; /* The x value at the current mesh node. */
    int n = p->nvars; /* The total number of dependant variables */

    matrix *InitSoln;
    InitSoln = CreateMatrix(p->nrows*p->nvars, 1);

    for(i=var; i<p->nrows*n; i+=n) {
        //x = valV(p->mesh->nodes, i/n);
        setval(InitSoln, value, i, 0);
    }

    return InitSoln;
}

Mesh1D* MoveMeshF(struct fe1d *p, Mesh1D *orig, double t,
                  double (*F)(struct fe1d *, double, double))
{
    int i;
    double x, dx, Fval;
    vector *new;

    new = CreateVector(len(orig->nodes));

    setvalV(new, 0, valV(orig->nodes, 0)); /* The first node doesn't change */

    for(i=1; i<len(orig->nodes); i++) {
        dx = valV(orig->nodes, i) - valV(orig->nodes, i-1);
        x = valV(orig->nodes, i);
        Fval = F(p, x, t);
        setvalV(new, i, Fval*dx+valV(new, i-1));
    }

    return Remesh1D(orig, new);
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

/**
 * Evaluate a one-dimensional solution at a particular local coordinate within
 * a specific element. This interpolates the already calculated solution using
 * the finite element basis functions. For problems with multiple dependent
 * variables, "var" specifies which one to interpolate. If the value supplied is
 * equal to -1, then the global x coordinate for that value of xi is returned
 * instead.
 * @param p Finite element problem structure
 * @param var Number corresponding to the dependent variable to interpolate
 * @param elem Mesh element to use when interpolating
 * @param s Set of dependent variable values at each node
 * @param xi Local coordinate value. Must be between 0 and 1.
 *
 * @returns Variable value at xi in the desired element.
 */
/* This function needs to be modified to support hermite cubic basis functions.
 */
double EvalSoln1D(struct fe1d *p, int var, Elem1D *elem, solution *s, double xi)
{   
    int i;
    double result = 0;
    int n = p->b->n;
    int nvars = p->nvars;

    /* Return the x (global) coordinate the corresponds to the xi (local)
     * coordinate when -1 is supplied for "var". */
    if(var == -1)
        for(i=0; i<n; i++)
            result += p->b->phi[i](xi) * valV(elem->points, i);

    /* Find the value of the desired variable at xi */
    for(i=0; i<n; i++) {
        result += p->b->phi[i](xi)
                  * val(s->val, valV(elem->map, i)*nvars+var, 0);
    }
    
    return result;
}

/**
 * Evaluate a one-dimensional solution at a given value of the global
 * independent variable, x.
 * @param p Finite element problem structure
 * @param var Number of the variable to interpolate
 * @param s Solution values
 * @param x Global x coordinate
 * @returns Interpolated solution at x
 * @see EvalSoln1D
 */
double EvalSoln1DG(struct fe1d *p, int var, solution *s, double x)
{
    int i;
    double x1, x2, xi, F, Fp, dx, h;
    Elem1D *e;
    e = NULL;
    /* Figure out which element the desired x value is in */
    /* If x is equal to the left endpoint of the mesh, return the first
     * element. */
    if(x == p->mesh->x1)
        e = p->mesh->elem[0];
    /* If it's equal to the right endpoint, return the last element. */
    else if(x == p->mesh->x2)
        e = p->mesh->elem[p->mesh->nelem-1];
    /* Otherwise, go through each element of the remaining elements one by one
     * until a match is found. Doing the previous two steps is likely
     * retundant. */
    else
            for(i=0; i<p->mesh->nelem; i++) {
                x1 = valV(p->mesh->elem[i]->points, 0);
                x2 = valV(p->mesh->elem[i]->points,
                                len(p->mesh->elem[i]->points));

                if(x >= x1 && x <= x2)
                    e = p->mesh->elem[i];
            }
    /* If we haven't found the element, quit the program. Something is likely
     * very wrong with the code. */
    if(!e) {
        printf("Failure to locate the element for point x = %g.\n", x);
        printf("Exiting.\n");
        exit(0);
    }
    
    /* Find the local coordinate, given the global coordinate using Newton's
     * method. */
    xi = .5; /* Since the range for xi is 0 to 1, just use 0.5 as the initial
              * guess */
    h = 1e-5; /* Tolerance for taking derivatives and Newton's method
               * convergence */
    do {
        F = EvalSoln1D(p, -1, e, s, xi);
        Fp= (EvalSoln1D(p, -1, e, s, xi+h)-EvalSoln1D(p, -1, e, s, xi-h))/(2*h);
        dx = -F/Fp;
        xi = xi + dx;
    } while(fabs(dx) > h);
    
    /* Return the desired value. */
    return EvalSoln1D(p, var, e, s, xi);
}

/**
 * Print out the nodal values at a given time step.
 * @param p Finite element problem structure
 * @param t Time step number
 */
void PrintSolution(struct fe1d *p, int t)
{
    solution *s;
    s = FetchSolution(p, t);
    mtxprnt(s->val);
}

/**
 * @brief Calculate the time derivative of the solution to a problem
 * @param problem The FE problem to use
 * @param x The calculated solution
 * @return The time derivative of the solution
 */
matrix* CalcTimeDerivative(struct fe1d *problem, matrix *x)
{
    matrix *tmp1, *tmp2, *dxdt;
    tmp1 = mtxmul(problem->J, x);
    mtxneg(tmp1);
    tmp2 = mtxadd(problem->F, tmp1);

    dxdt = SolveMatrixEquation(problem->dJ, tmp2);

    DestroyMatrix(tmp1);
    DestroyMatrix(tmp2);

    return dxdt;
}


