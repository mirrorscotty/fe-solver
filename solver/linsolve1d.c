/**
 * @file linsolve1d.c
 * Several linear finite element solvers for linear, one-dimensional problems.
 */

#include <math.h>
#include <stdio.h>
#include "matrix.h"
#include "finite-element1d.h"
#include "solve.h"

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
 * Implicit time integration solver. (Linear version)
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
    matrix *u, /* Used to store the solution at the previous time step */
           *du; /* Derivative of u at the previous time step */
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

    /* Apply boundary conditions here so that any Dirchlet boundary conditions
     * are imposed properly for the transient solver. Doing it here, rather than
     * above guarantees that the specified value on the boundaries doesn't
     * change randomly. */
    DestroyMatrix(problem->J);
    DestroyMatrix(problem->F);
    problem->J = tmp2;
    problem->F = tmp3;
    problem->applybcs(problem);

    /* Solve the matrix equation to find the nodal values */
    result = SolveMatrixEquation(problem->J, problem->F);

    /* Solve for the time derivative of the result. */
    dresult = CalcTimeDerivative(problem, result);

    StoreSolution(problem, result, dresult);
    /* Delete the time derivative of the pervious solution to save memory. */
    DeleteTimeDeriv(prev);

    return result;
}

