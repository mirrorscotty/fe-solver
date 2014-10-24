/**
 * @file nlinsolve2d.c
 * Nonlinear, 2D finite element solvers for static problems only.
 */

#include <math.h>
#include <stdio.h>
#include "matrix.h"
#include "finite-element.h"
#include "solve.h"

/**
 * @brief Check the convergance of the nonlinear solver.
 *
 * @param problem A pointer to the problem struct
 * @param dt A matrix containing the change in the x vector
 * @returns 0 if not converged, and 1 if we're done iterating.
 */
int CheckConverg(struct fe *problem, matrix *dx)
{
    int rows = nRows(dx);
    int i;
    
    for(i=0; i<rows; i++) {
        if(fabs(val(dx, i, 0)) > problem->tol)
            return 0;
    }
    return 1;
}

/**
 * @brief Nonlinear finite element solver.
 *
 * This does all the same stuff as the linear
 * solver (assembling the global matricies), but utilizes an initial guess to
 * calculate them. It then iterates using Newton's method and updates this
 * guess at each iteration.
 *
 * @param problem The struct containing the problem to solve
 * @param guess The initial guess
 *
 * @returns A matrix with the solution at each node
 */
matrix* NLinSolve(struct fe *problem, matrix *guess)
{
    matrix *dx; /* How much to update the guess by */
    matrix *newguess;
    //int rows = (2*problem->mesh->nelemx+1)*(2*problem->mesh->nelemy+1);
    int rows = problem->nrows;
    int iter = 0;
    int maxiter = 500;
    
    if(!guess) {
        guess = CreateMatrix(rows*problem->nvars, 1);
    }

    do {
        iter++;
        
        if(problem->J)
            DestroyMatrix(problem->J);
        if(problem->F)
            DestroyMatrix(problem->F);
        
        AssembleJ(problem, guess);
        problem->F = CreateMatrix(rows*problem->nvars, 1);
        //AssembleF(problem, guess);
        problem->applybcs(problem);

        //if(!CalcDeterminant(problem->J)) {
        //    iter = -1;
        //   break;
        //}
        
        CalcResidual(problem, guess);
        dx = SolveMatrixEquation(problem->J, problem->R);
        newguess = mtxadd(guess, dx);
        DestroyMatrix(guess);
        guess = newguess;

        /* Quit if we've reached the maximum number of iterations */
        if(iter == maxiter)
            break;
        
        printf("\rIteration %d", iter); // Print the current iteration number to the console.
        fflush(stdout); // Flush the output buffer.
        
    } while(!CheckConverg(problem, dx));
    /* ^^ Also quit if the dx variable is small enough. */
    
    if(iter == -1)
        /* If we've determined the matrix to be singular by calculating the
         * determinant, then output the appropriate error message. */
        printf("\rSingular matrix.\n");
    else if(iter == maxiter)
        /* If the solver didn't find a solution in the specified number of
         * iterations, then say so. */
        printf("\rNonlinear solver failed to converge. Maximum number of iterations reached.\n");
    else
        /* Print out the number of iterations it took to converge
         * successfully. */
        printf("\rNonlinear solver converged after %d iterations.\n", iter);
    
    return guess;
}

