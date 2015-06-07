/**
 * @file auxsoln.c
 * Functions to allow for the solving of ODEs based on the solution to a PDE
 * that was calculated using finite element method.
 */

#include <stdlib.h>

#include "auxsoln.h"
#include "finite-element1d.h"
#include "matrix.h"

/**
 * Set up the finite element problem structure to accomodate solving ODEs based
 * on the PDE solutions.
 * @param p Finite element problem structure
 * @param nvars Number of variables (equations) to solve.
 */
void FE1DInitAuxSolns(struct fe1d *p, int nvars)
{
    int i;
    p->extravars = nvars;
    p->auxsolns = (solution***) calloc(nvars, sizeof(solution**));
    for(i=0; i<nvars; i++)
        p->auxsolns[i] = (solution**) calloc(p->maxsteps, sizeof(solution*));
    return;
}

/**
 * Integrate a time-dependent ODE whose solution is based on the PDE solution.
 * This function uses Euler's method to perform time integration.
 * @param p Finite element problem structure
 * @param var Number of the variable to use as an input for the ODE
 * @param extravar Output variable for the ODE solution
 * @param eqn Equation solved for DuDt, where u is the dependent variable of
 *       the ODE.
 */
void SolveODE(struct fe1d *p, int var, int extravar,
              double (*eqn)(double, double, double), double IC)
{
    int t, i;
    int tmax = p->maxsteps;
    solution *s, *tmp, *tmpprev;

    double dt, T, cprev;
    int rows;

    s = FetchSolution(p, 0);
    rows = nRows(s->val); /* Determine how many nodes there are */

    /* Create the initial condition */
    tmp = CreateSolution(0, 0, CreateMatrix(rows, 1));
    for(i=0; i<rows; i++)
        setval(tmp->val, IC, i, 0);
    p->auxsolns[extravar][0] = tmp;

    for(t=1; t<tmax; t++) {
        s = FetchSolution(p, t);
        dt = uscaleTime(p->charvals, s->dt);

        tmpprev = tmp;
        tmp = CreateSolution(t, dt, CreateMatrix(rows, 1));
        for(i=0; i<rows; i++) {
            cprev = val(tmpprev->val, i, 0);
            T = uscaleTemp(p->charvals, val(s->val, i, 0));

            setval(tmp->val, eqn(cprev, T, dt), i, 0);
        }
        p->auxsolns[extravar][t] = tmp;
    }

    return;
}

/**
 * Print out an ODE solution.
 * @param p Finite element problem structure
 * @param var Number of the ODE to print the solution for.
 * @param t Print out the solution at this time step for all values of x.
 */
void PrintAuxSoln(struct fe1d *p, int var, int t)
{
    solution *s;
    s = p->auxsolns[var][t];
    mtxprnt(s->val);
}

/**
 * Return an ODE solution structure
 * @param p Finite element problem structure
 * @param var Number of the ODE to return the solution for
 * @param t Return the solution at all values of x for time t
 */
solution* FetchAuxSoln(struct fe1d *p, int var, int t)
{
    return p->auxsolns[var][t];
}

void InitAuxStep(struct fe1d *p, int var, int t)
{
    int nnodes = len(p->mesh->nodes);
    matrix *m;
    solution *s;

    m = CreateMatrix(nnodes, 1);
    s = CreateSolution(t, -1, m);

    p->auxsolns[var][t] = s;
}

