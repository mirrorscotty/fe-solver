/**
 * @file linsolve2d.c
 * Two-dimensional finite element solver for linear, static problems
 */
#include <math.h>
#include <stdio.h>
#include "matrix.h"
#include "finite-element.h"
#include "solve.h"

/**
 * @brief Linear finite element solver.
 * It assembles the global jacobian and load vectors, then solves the whole
 * system. It returns a matrix with the solution.
 * @param problem The struct representing the problem to solve
 * @returns A matrix with the solution at each node
 */
matrix* LinSolve(struct fe *problem)
{
    matrix *guess;
    int rows = problem->nrows;
    //int rows = (problem->mesh->nelemx+1)*(problem->mesh->nelemy+1);
    guess = CreateMatrix(rows*problem->nvars, 1);
    AssembleJ(problem, guess);
    problem->F = CreateMatrix(rows*problem->nvars, 1);
    problem->applybcs(problem);

    /* Calculate the determinate of the jacobian matrix to see if we're going to
     * good results. Currently this takes forever, and is commented for that
     * reason. */
//    if(fabs(CalcDeterminant(problem->J)) - 1e-10 < 0)
//        puts("Singular Matrix.");

    return SolveMatrixEquation(problem->J, problem->F);
}

