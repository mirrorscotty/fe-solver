#include <math.h>
#include <stdio.h>
#include "matrix.h"
#include "finite-element.h"
#include "finite-element1d.h"
#include "mtxsolver.h"
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

int CheckConverg1D(struct fe1d *problem, matrix *dx)
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

/* Todo: replace this function and the above one by a macro for simplicity. */
/**
 * @brief 1D version of LinSolve.
 * @see LinSolve
 */
matrix* LinSolve1D(struct fe1d *problem)
{
    matrix *guess;
    int rows = problem->nrows;
    guess = CreateMatrix(rows*problem->nvars, 1);
    AssembleJ1D(problem, guess);
    problem->F = CreateMatrix(rows*problem->nvars, 1);
    problem->applybcs(problem);
    
    /* Calculate the determinate of the jacobian matrix to see if we're going to
     * good results. Currently this takes forever, and is commented for that
     * reason. */
//    if(fabs(CalcDeterminant(problem->J)) - 1e-10 < 0)
//        puts("Singular Matrix.");

    return SolveMatrixEquation(problem->J, problem->F);
}

/**
 * @brief Explicit time integration algorithm (Forward Euler)
 * 
 * This solver fails horribly if a Neumann boundary condition is imposed. Also,
 * it may not be entirely stable anyway.
 *
 * @param problem The struct containing the problem to solve
 * @returns A matrix with the solutions at each node
 */
matrix* LinSolve1DTrans(struct fe1d *problem)
{
    matrix *du, *u;
    matrix *tmp1, *tmp2, *tmp3, *tmp4;
    matrix *dresult, *result;
    solution *prev; /* The solution at the previous time step */

    /* Initialize the matricies */
    DestroyMatrix(problem->J);
    DestroyMatrix(problem->dJ);
    DestroyMatrix(problem->F);
    AssembleJ1D(problem, NULL);
    AssembledJ1D(problem, NULL);
    AssembleF1D(problem, NULL);
    problem->applybcs(problem);


    /* Get the previous solution */
    prev = FetchSolution(problem, problem->t-1);
    du = prev->dval;
    u = prev->val;

    tmp1 = mtxmulconst(du, problem->dt);
    tmp2 = mtxadd(u, tmp1);
    mtxneg(tmp2);
    tmp3 = mtxmul(problem->J, tmp2);

    tmp4 = mtxadd(problem->F, tmp3);

    /* Solve for du/dt */
    dresult = SolveMatrixEquation(problem->dJ, tmp4);

    DestroyMatrix(tmp1);
    DestroyMatrix(tmp2);
    DestroyMatrix(tmp3);
    DestroyMatrix(tmp4);

    /* Solve for u */
    tmp1 = mtxmul(problem->dJ, dresult);
    mtxneg(tmp1);
    tmp2 = mtxadd(problem->F, tmp1);
    result = SolveMatrixEquation(problem->J, tmp2);

    DestroyMatrix(tmp1);
    DestroyMatrix(tmp2);

    StoreSolution(problem, result, dresult);
    
    return result;
}

/**
 * @brief Implicit time integration solver.
 *
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
 * @param problem The problem to solve
 * @returns A matrix with the solutions at each node
 */
matrix* LinSolve1DTransImp(struct fe1d *problem)
{
    /* Constants that determine the integration algorithm. If a=1 and b=0, then
     * the backward difference method is being used. For the trapazoid rule,
     * a=2, and b=-1. */
    int a, b;
    matrix *du, *u;
    /* Used to store intermediate calculations */
    matrix *tmp1, *tmp2, *tmp3, *tmp4;
    matrix *dresult, *result;
    solution *prev; /* Solution at previous time step */

    /* Use BD for testing */
    a = 1; b = 0;

    /* Initialize the matricies */
    DestroyMatrix(problem->J);
    DestroyMatrix(problem->dJ);
    DestroyMatrix(problem->F);
    AssembleJ1D(problem, NULL);
    AssembledJ1D(problem, NULL);
    AssembleF1D(problem, NULL);
    problem->applybcs(problem);

    /* Get the previous solution */
    prev = FetchSolution(problem, problem->t-1);
    du = prev->dval;
    u = prev->val;

    /* Calculate the matrix to be multiplied by the unknown vector */
    tmp1 = mtxmulconst(problem->dJ, a/problem->dt);
    tmp2 = mtxadd(tmp1, problem->J);

    DestroyMatrix(tmp1);
    
    /* Determine the right hand side of the equation. The if statement is there
     * so that this function can be used for both TR and BD, but not waste an
     * extra step multiplying by 0 for BD. */
    tmp1 = mtxmulconst(problem->dJ, a/problem->dt);
    tmp3 = mtxmul(tmp1, u);

    DestroyMatrix(tmp1);

    if(b) {
        tmp1 = mtxmulconst(problem->dJ, b);
        tmp4 = mtxmul(tmp1, du);
        mtxneg(tmp4);
        DestroyMatrix(tmp1);

        tmp1 = mtxadd(tmp3, tmp4);
        DestroyMatrix(tmp4);
        DestroyMatrix(tmp3);
    } else {
        tmp1 = tmp3;
    }

    tmp3 = mtxadd(tmp1, problem->F);
    DestroyMatrix(tmp1);

    result = SolveMatrixEquation(tmp2, tmp3);

    DestroyMatrix(tmp2);
    DestroyMatrix(tmp3);


    /* Solve for the time derivative of the result. */
    dresult = CalcTimeDerivative(problem, result);

    StoreSolution(problem, result, dresult);
    /* Delete the time derivative of the pervious solution to save memory. */
    DeleteTimeDeriv(prev);

    if(problem->t == 2)
    mtxprnt(result);

    return result;
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

/**
 * @brief Nonlinear finite element solver. (1D version)
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

        printf("\rIteration %d", iter); // Print the current iteration number to the console.
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
        printf("\rNonlinear solver failed to converge. Maximum number of iterations reached.\n");
    else
        /* Print out the number of iterations it took to converge
         * successfully. */
        printf("\rNonlinear solver converged after %d iterations.\n", iter);

    return guess;
}

matrix* NLinSolve1DTransImp(struct fe1d *problem, matrix *guess)
{
    int iter = 0;
    int maxiter = 50;

    int a, b;
    solution *prev;
    matrix *J, *R;
    matrix *tmp1, *tmp2, *tmp3, *tmp4, *tmp5, *tmp6;
    matrix *u, *du, *dx;
    matrix *dguess; /* Time derivative of the guess */

    /* Use BD */
    a = 1; b = 0;

    /* Get the previous solution */
    prev = FetchSolution(problem, problem->t-1);
    du = prev->dval;
    u = prev->val;

    /* Set the initial guess to the solution from the previous time step if it
     * isn't supplied. */
    if(!guess)
        guess = CopyMatrix(u);

    do {
        iter++;

        /* Quit if we've reached the maximum number of iterations */
        if(iter == maxiter)
            break;

        problem->guess = guess;

        /* Initialize the matricies */
        DestroyMatrix(problem->J);
        DestroyMatrix(problem->dJ);
        DestroyMatrix(problem->F);
        AssembleJ1D(problem, guess);
        AssembledJ1D(problem, guess);
        AssembleF1D(problem, guess);
        problem->applybcs(problem);

        problem->guess = NULL;

        /* Calculate the Jacobian matrix */
        tmp1 = mtxmulconst(problem->dJ, a/problem->dt);
        J = mtxadd(tmp1, problem->J);

        DestroyMatrix(tmp1);

        /* Calculate the residual matrix */
        tmp1 = mtxmul(J, guess);
        mtxneg(tmp1);
        tmp2 = mtxadd(tmp1, problem->F);
        tmp3 = mtxmul(problem->dJ, u);
        tmp4 = mtxmulconst(tmp3, a/problem->dt);
        tmp5 = mtxadd(tmp2, tmp4);
        if(b) {
            tmp6 = mtxmulconst(du, b);
            R = mtxadd(tmp5, tmp6);
            DestroyMatrix(tmp6);
        } else {
            R = CopyMatrix(tmp5);
        }
        DestroyMatrix(tmp1);
        DestroyMatrix(tmp2);
        DestroyMatrix(tmp3);
        DestroyMatrix(tmp4);
        DestroyMatrix(tmp5);

        /* Solve for dx */
        //mtxneg(R);
        dx = SolveMatrixEquation(J, R);
        mtxprnt(guess);

        /* Add the change in the unknowns to the guess from the previous
         * iteration */
        tmp1 = mtxadd(guess, dx);
        DestroyMatrix(guess);
        guess = tmp1;

        printf("\rIteration %d", iter); // Print the current iteration number to the console.
        fflush(stdout); // Flush the output buffer.
    } while(!CheckConverg1D(problem, dx));

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

    /* Solve for the time derivative of the result. */
    dguess = CalcTimeDerivative(problem, guess);

    StoreSolution(problem, guess, dguess);
    /* Delete the time derivative of the pervious solution to save memory. */
    DeleteTimeDeriv(prev);

    return guess;
}
