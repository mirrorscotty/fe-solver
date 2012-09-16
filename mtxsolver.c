#include <math.h>
#include <stdio.h>
#include "matrix.h"
#include "finite-element.h"
#include "finite-element1d.h"

void SwapRows(matrix*, int, int);
int FindPivot(matrix*, int);

/* Code from http://compprog.wordpress.com/2007/12/11/gaussian-elimination/ */

void ForwardSubstitution(matrix* a) {
	int i, j, k, max;
	int n = nRows(a);
	double t;
	for (i = 0; i < n; ++i) {
		max = i;
		for (j = i + 1; j < n; ++j)
			if (fabs(val(a, j, i)) > fabs(val(a, max, i)))
				max = j;
		
		for (j = 0; j < n + 1; ++j) {
			t = val(a, max, j);
			setval(a, val(a, i, j), max, j);
			setval(a, t, i, j);
		}
		
		for (j = n; j >= i; --j)
			for (k = i + 1; k < n; ++k)
				setval(a, val(a, k, j) - val(a, k, i)/val(a, i, i) * val(a, i, j), k, j);

	}
}

void ReverseElimination(matrix *a) {
	int i, j;
	int n = nRows(a);
	for (i = n - 1; i >= 0; --i) {
		setval(a, val(a, i, n)/val(a, i, i), i, n);
		setval(a, 1, i, i);
		for (j = i - 1; j >= 0; --j) {
			setval(a, val(a, j, n) - val(a, j, i) * val(a, i, n), j, n);
			setval(a, 0, j, i);
		}
	}
}

/* Solve an matrix equation of the form Ax=B, where x is the vector of
 * unknowns. */
matrix* SolveMatrixEquation(matrix *A, matrix *B)
{
	matrix *C, *u;
	C = AugmentMatrix(A, B);
	ForwardSubstitution(C);
	ReverseElimination(C);
	u = ExtractColumn(C, nCols(A));
	DestroyMatrix(C);
	return u;
}

/* Check the convergance of the nonlinear solver. Return 0 if not converged,
 * and 1 if we're done iterating. */
int CheckConverg(struct fe *problem, matrix *dx)
{
    int rows = nRows(dx);
    int i;
    
    for(i=0; i<rows; i++) {
        if(val(dx, i, 0) > problem->tol)
            return 0;
    }
    return 1;
}
    
/* Linear finite element solver. It assembles the global jacobian and load
 * vectors, then solves the whole system. It returns a matrix with the
 * solution. */
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

/* 1D version of the above function. */
/* Todo: replace this function and the above one by a macro for simplicity. */
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

/* Explicit time integration algorithm (Forward Euler)
 * This solver fails horribly if a Neumann boundary condition is imposed. Also,
 * it may not be entirely stable anyway. */
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

/* Implicit time integration solver. Works with both backward difference
 * integration (BD) and the trapazoid rule (TR). In either case, the algorithm
 * is unconditionally stable. */
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
    tmp1 = mtxmul(problem->J, result);
    mtxneg(tmp1);
    tmp2 = mtxadd(problem->F, tmp1);

    dresult = SolveMatrixEquation(problem->dJ, tmp2);

    DestroyMatrix(tmp1);
    DestroyMatrix(tmp2);

    StoreSolution(problem, result, dresult);

    return result;
}

/* Nonlinear finite element solver. This does all the same stuff as the linear
 * solver (assembling the global matricies), but utilizes an initial guess to
 * calculate them. It then iterates using Newton's method and updates this
 * guess at each iteration. */
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
        
        if(iter == maxiter)
            break;
        
        printf("\rIteration %d", iter); // Print the current iteration number to the console.
        fflush(stdout); // Flush the output buffer.
        
    } while(!CheckConverg(problem, dx));
    
    if(iter == -1)
        printf("\rSingular matrix.\n");
    else if(iter == maxiter)
        printf("\rNonlinear solver failed to converge. Maximum number of iterations reached.\n");
    else
        printf("\rNonlinear solver converged after %d iterations.\n", iter);
    
    return guess;
}
