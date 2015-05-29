/**
 * @file nlinsolve1d.c
 * Finite element solvers for both one-dimensional problems. Includes solvers
 * for static and transient problems.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"
#include "finite-element1d.h"
#include "solve.h"

/**
 * Check the convergence of the nonlinear solver for a 1D problem.
 *
 * @param problem A pointer to the problem struct
 * @param dt A matrix containing the change in the x vector
 * @returns 0 if not converged, and 1 if we're done iterating.
 *
 * @see CheckConverg
 */
int CheckConverg1D(struct fe1d *problem, matrix *dx)
{
    int rows = nRows(dx);
    int i;

    for(i=0; i<rows; i++) {
        /* Chech to ensure that each element of the dx matrix is less than the
         * tolerance. If not, then return 0, indicating more iterations are
         * required. */
        if(fabs(val(dx, i, 0)) > problem->tol)
            return 0;
        
        /* If one of the values in the dx matrix is NaN, then clearly something
         * is wrong. The best thing to do is to quit and spit out an error. */
        if(isnan(val(dx, i, 0))) {
            printf("Nonlinear solver failed to converge."
                       "Solver returned value of \"NaN\".\n"
                       "Failed to calculate solution at time step %d of %d\n"
                       "Exiting.\n",
                   problem->t, problem->maxsteps);
            exit(0);
        }
    }
    
    /* All the values are less than the specified tolerance, so we're good to
     * go. */
    return 1;
}

/**
 * Nonlinear finite element solver. (1D version)
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
matrix* NLinSolve1D(struct fe1d *problem, matrix *guess)
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

        AssembleJ1D(problem, guess);
        problem->F = CreateMatrix(rows*problem->nvars, 1);
        //AssembleF(problem, guess);
        problem->applybcs(problem);

        //if(!CalcDeterminant(problem->J)) {
        //    iter = -1;
        //   break;
        //}

        /* ToDo: Write this function! */
        CalcResidual1D(problem, guess);
        dx = SolveMatrixEquation(problem->J, problem->R);
        newguess = mtxadd(guess, dx);
        DestroyMatrix(guess);
        guess = newguess;

        /* Quit if we've reached the maximum number of iterations */
        if(iter == maxiter)
            break;

        /* Print the current iteration number to the console. */
        printf("\rIteration %d", iter); 
        fflush(stdout); // Flush the output buffer.

    } while(!CheckConverg1D(problem, dx));
    /* ^^ Also quit if the dx variable is small enough. */

    if(iter == -1)
        /* If we've determined the matrix to be singular by calculating the
         * determinant, then output the appropriate error message. */
        printf("\rSingular matrix.\n");
    else if(iter == maxiter)
        /* If the solver didn't find a solution in the specified number of
         * iterations, then say so. */
        printf("\rNonlinear solver failed to converge."
               "Maximum number of iterations reached.\n");
    else
        /* Print out the number of iterations it took to converge
         * successfully. */
        printf("\rNonlinear solver converged after %d iterations.\n", iter);

    return guess;
}

/**
 * Implicit time integration solver. (Nonlinear version)
 * Works with both backward difference
 * integration (BD) and the trapazoid rule (TR). In either case, the algorithm
 * is supposed to be unconditionally stable. This function only solves the
 * problem at the next time step. In order to solve for the desired variables
 * over the entire time domain, call this function repeatedly.
 *
 * The time deriviative is calculated using the following formula:
 * y'(t1) = a/dt * [y(t1) - y(t0)] + b*y'(t0)
 * Here, t1 is the time at the current time step and t0 is at the previous one.
 * The variables a and b are constants that determine which algorithm is used
 * to approximate the derivative. For a=1, b=0, backward difference is used,
 * and for a=2, b=-1, trapazoid rule is used.
 *
 * @param problem Finite element problem to solve
 * @param guess Initial guess at the solution for the next time step. If this
 *      is not supplied (NULL), then the solution is predicted based on the
 *      previous time steps.
 * @returns Pointer to a matrix of the calculated values. This can be safely
 *      disregarded since the matrix is also stored in the finite element
 *      problem structure.
 */
matrix* NLinSolve1DTransImp(struct fe1d *problem, matrix *guess)
{
    int iter = 0; /* Current iteration number */
    int maxiter = 5000; /* Maximum allowed iterations before terminating */

    /* Constants that determine the integration algorithm. If a=1 and b=0, then
     * the backward difference method is being used. For the trapazoid rule,
     * a=2, and b=-1. */
    int a, b;
    solution *prev; /* Solution at the previous time step */
    matrix *J, *F, *R; /* Jacobian matrix and the residuals matrix */
    matrix *tmp1, *tmp2; /* Temporary matricies */
    matrix *u, *du, *dx;
    matrix *dguess; /* Time derivative of the guess */

    /* Use BD */
    a = 1; b = 0;
    /* Note: This solver doesn't work at all with anything but backward
     * difference for time integration. The linear one doesn't either, but it
     * at least contains some code for it. */

    /* Get the previous solution */
    prev = FetchSolution(problem, problem->t-1);
    du = prev->dval;
    u = prev->val;

    /* Predict the next solution if an initial guess isn't supplied. */
    if(!guess)
        guess = PredictSolnO0(problem);

    dx = NULL;
    //exit(0);

    do {
        iter++;

        /* Quit if we've reached the maximum number of iterations */
        if(iter == maxiter)
            break;

        problem->guess = guess;

        /* Delete the matricies from the previous iteration */
        DestroyMatrix(problem->J);
        DestroyMatrix(problem->dJ);
        DestroyMatrix(problem->F);
        /* Delete the dx matrix from the previous iteration. */
        if(dx)
            DestroyMatrix(dx); 
        /* Initialize the matricies */
        AssembleJ1D(problem, guess);
        AssembledJ1D(problem, guess);
        AssembleF1D(problem, guess);

        //problem->guess = NULL;

        /* Calculate the Jacobian matrix */
        tmp1 = mtxmulconst(problem->dJ, a/problem->dt);
        J = mtxadd(tmp1, problem->J);
        DestroyMatrix(tmp1);

        /* Calculate the load vector */
        tmp1 = mtxmulconst(problem->dJ, a/problem->dt);
        tmp2 = mtxmul(tmp1, u);
        F = mtxadd(problem->F, tmp2);
        DestroyMatrix(tmp1);
        DestroyMatrix(tmp2);

        /* Apply boundary conditions. */
        DestroyMatrix(problem->J);
        DestroyMatrix(problem->F);
        problem->J = J;
        problem->F = F;
        problem->applybcs(problem);

        /* Calculate the residual vector */
        tmp1 = mtxmul(problem->J, guess);
        mtxneg(tmp1);
        R = mtxadd(problem->F, tmp1);
        DestroyMatrix(tmp1);

        //mtxprnt(problem->J);
        //puts("");
        //mtxprnt(problem->dJ);
        //puts("");
        //mtxprnt(problem->F);

        /* Solve for dx */
        dx = SolveMatrixEquation(J, R);

        /* Delete the residual matrix and the Jacobian matrix*/
        DestroyMatrix(R);

        /* Add the change in the unknowns to the guess from the previous
         * iteration */
        tmp1 = mtxadd(guess, dx);
        DestroyMatrix(guess);
        guess = tmp1;

#ifdef VERBOSE_OUTPUT
        /* Print the current iteration number to the console. */
        printf("\rIteration %d", iter);
        fflush(stdout); // Flush the output buffer.
#endif
    } while(!CheckConverg1D(problem, dx));
    /* Delete the final dx matrix */
    DestroyMatrix(dx);

    if(iter==maxiter) {
        printf("Nonlinear solver failed to converge. Maximum number of iterations reached.\nFailed to calculate solution at time step %d of %d\nExiting.\n",
                problem->t, problem->maxsteps);
        exit(0);
    }

#ifdef VERBOSE_OUTPUT
    if(iter == -1)
        /* If we've determined the matrix to be singular by calculating the
         * determinant, then output the appropriate error message. */
        printf("\rSingular matrix.\n");
    else if(iter == maxiter)
        /* If the solver didn't find a solution in the specified number of
         * iterations, then say so. */
        printf("\rNonlinear solver failed to converge."
               "Maximum number of iterations reached.\n");
    else
        /* Print out the number of iterations it took to converge
         * successfully. */
        printf("\rNonlinear solver converged after %d iterations.\n", iter);
#endif

    /* Solve for the time derivative of the result. */
    dguess = CalcTimeDerivative(problem, guess);

    StoreSolution(problem, guess, dguess);
    /* Delete the time derivative of the pervious solution to save memory. */
    //DeleteTimeDeriv(prev);

    return guess;
}

