/**
 * @file finite-element1d.h
 */

#ifndef FINITE_ELEMENT1D_H
#define FINITE_ELEMENT1D_H

#include "matrix.h"
#include "basis.h"
#include "mesh1d.h"
#include "solution.h"

#include "scaling_ht.h"

struct fe1d;

/**
 * Structure used to define a 1-dimensional finite element problem.
 */
struct fe1d {
    basis *b; /**< Set of interpolation functions to use when solving */
    Mesh1D *mesh; /**< Mesh that is used to discretize the domain */

    matrix *J; /**< Jacobian matrix*/
    matrix *dJ; /**< Portion of Jacobian multiplying time-dependent part of
                    the residual */
    matrix *R; /**< Residual matrix */
    matrix *F; /**< Load vector */

    matrix *guess;

    solution **soln; /**< Stored solutions for each time step. */
    int t; /**< Current time index. */
    int maxsteps; /**< Maximum number of time steps to be allowed. */
    double dt; /**< The time step for the initial step. */

    int nvars; /**< Number of dependant variables. */
    int nrows; /**< Number of rows (nodes) in the final solution. */

    double tol; /**< Tolerance for the nonlinear solver */

    /* Functions to make the element mass/stiffness matrix and the load
     * vector */
    matrix* (*makedj)(struct fe1d*, Elem1D*, matrix*);
    matrix* (*makej)(struct fe1d*, Elem1D*, matrix*);
    matrix* (*makef)(struct fe1d*, Elem1D*, matrix*);

    void (*applybcs)(struct fe1d*); /**< Function to apply all the relevant
                                     *   boundary conditions for the problem. */

    /* Stuff for ODEs */
    solution ***auxsolns;
    int extravars; /**< Number of ODE variables to solve for */

    /* Scaling stuff */
    scaling_ht charvals; /**< Scaling values for heat transfer. This structure
                          *   holds the characteristic length, temperature,
                          *   thermal diffusivity, etc. used to make the
                          *   heat transfer portion of the problem
                          *   dimensionless.
                          */
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

Mesh1D* MoveMeshF(struct fe1d *, Mesh1D *, double,
        double (*F)(struct fe1d *, double, double));

int StoreSolution(struct fe1d*, matrix*, matrix*);
solution* FetchSolution(struct fe1d*, int);
double EvalSoln1D(struct fe1d*, int, Elem1D*, solution*, double);
void PrintSolution(struct fe1d*, int);

matrix* CalcTimeDerivative(struct fe1d*, matrix*);

#endif

