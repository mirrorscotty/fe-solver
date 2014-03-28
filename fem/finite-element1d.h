#ifndef FINITE_ELEMENT1D_H
#define FINITE_ELEMENT1D_H

#include "matrix.h"
#include "basis.h"
#include "mesh1d.h"
#include "solution.h"

#include "scaling_ht.h"

struct fe1d;

struct fe1d {
    basis *b; /* Set of interpolation functions to use when solving */
    Mesh1D *mesh; /* Mesh that discretizes the domain */

    matrix *J; /* Jacobian */
    matrix *dJ; /* Portion of Jacobian multiplying time-dependent part of
                   the residual */
    matrix *R; /* Residual */
    matrix *F; /* Load vector */

    matrix *guess;

    solution **soln; /* Stored solutions for each time step. */
    int t; /* Current time index. */
    int maxsteps; /* Maximum number of time steps to be allowed. */
    double dt; /* The time step for the initial step. */

    int nvars; /* Number of dependant variables. */
    int nrows; /* Number of rows (nodes) in the final solution. */

    double tol; /* Tolerance for the nonlinear solver */

    /* Functions to make the element mass/stiffness matrix and the load vector */
    matrix* (*makedj)(struct fe1d*, Elem1D*, matrix*);
    matrix* (*makej)(struct fe1d*, Elem1D*, matrix*);
    matrix* (*makef)(struct fe1d*, Elem1D*, matrix*);

    /* Function to apply all the relevant boundary conditions for the problem. */
    void (*applybcs)(struct fe1d*);

    /* Stuff for ODEs */
    solution ***auxsolns;
    int extravars;

    /* Scaling stuff */
    scaling_ht charvals;
    scaling_ht chardiff;
};

struct fe1d* CreateFE1D(basis*, Mesh1D*,
                        matrix* (*)(struct fe1d*, Elem1D*, matrix*),
                        matrix* (*)(struct fe1d*, Elem1D*, matrix*),
                        matrix* (*)(struct fe1d*, Elem1D*, matrix*),
                        void (*)(struct fe1d*), int);
void DestroyFE1D(struct fe1d*);

matrix* AssembleJ1D(struct fe1d*, matrix*);
matrix* AssembledJ1D(struct fe1d*, matrix*);
matrix* AssembleF1D(struct fe1d*, matrix*);

matrix* CalcResidual1D(struct fe1d*, matrix*);

void FE1DTransInit(struct fe1d*, matrix*);
matrix* AssembleF1DTrans(struct fe1d*, matrix*);

void ApplyEssentialBC1D(struct fe1d*, int, int (*)(struct fe1d*, int),
                        double (*)(struct fe1d*, int));
void ApplyNaturalBC1D(struct fe1d*, int, int (*)(struct fe1d*, int),
                      double (*)(struct fe1d*, int));

matrix* GenerateInitialCondition(struct fe1d*, int, double (*)(double));
matrix* GenerateInitCondConst(struct fe1d*, int, double);

int StoreSolution(struct fe1d*, matrix*, matrix*);
solution* FetchSolution(struct fe1d*, int);
void PrintSolution(struct fe1d*, int);

#endif

