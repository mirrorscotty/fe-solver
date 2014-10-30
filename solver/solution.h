/**
 * @file solution.h
 * Storage for the calculated solutions to the PDEs at each time step.
 */
#ifndef SOLUTION_H
#define SOLUTION_H

#include "matrix.h"

/**
 * Structure to store the values solved for using finite element method.
 */
typedef struct {
    int t; /**< Time step number */
    double dt; /**< Time step size (Currently this shouldn't change from
                *   solution to solution if you want to avoid problems) */
    matrix *val; /**< The values of the variable(s) at the indicated time. */
    matrix *dval; /**< The values of the first time derivative. */
} solution;

solution* CreateSolution(int, double, matrix*);
void DestroySolution(solution*);
void DeleteTimeDeriv(solution*);

#endif

