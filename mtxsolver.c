#include <math.h>
#include <stdio.h>
#include "matrix.h"
#include "finite-element.h"

void SwapRows(matrix*, int, int);
int FindPivot(matrix*, int);

/* Code from http://compprog.wordpress.com/2007/12/11/gaussian-elimination/ */

void ForwardSubstitution(matrix* a) {
	int i, j, k, max;
	int n = mtxlen2(a);
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
	int n = mtxlen2(a);
	for (i = n - 1; i >= 0; --i) {
		setval(a, val(a, i, n)/val(a, i, i), i, n);
		setval(a, 1, i, i);
		for (j = i - 1; j >= 0; --j) {
			setval(a, val(a, j, n) - val(a, j, i) * val(a, i, n), j, n);
			setval(a, 0, j, i);
		}
	}
}

matrix* SolveMatrixEquation(matrix *A, matrix *B)
{
	matrix *C, *u;
	C = AugmentMatrix(A, B);
	ForwardSubstitution(C);
	ReverseElimination(C);
	u = ExtractColumn(C, mtxlen1(A));
	DestroyMatrix(C);
	return u;
}

/* Return 0 if not converged, and 1 if we're done iterating. */
int CheckConverg(struct fe *problem, matrix *dx)
{
    int rows = mtxlen2(dx);
    int i;
    
    for(i=0; i<rows; i++) {
        if(val(dx, i, 0) > problem->tol)
            return 0;
    }
    return 1;
}
    
matrix* LinSolve(struct fe *problem)
{
    matrix *guess;
    int rows = (2*problem->mesh->nelemx+1)*(2*problem->mesh->nelemy+1);
    //int rows = (problem->mesh->nelemx+1)*(problem->mesh->nelemy+1);
    guess = CreateMatrix(rows*problem->nvars, 1);
    AssembleJ(problem, guess);
    problem->F = CreateMatrix(rows*problem->nvars, 1);
    problem->applybcs(problem);
    
//    if(fabs(CalcDeterminant(problem->J)) - 1e-10 < 0)
//        puts("Singular Matrix.");

    return SolveMatrixEquation(problem->J, problem->F);
}

matrix* NLinSolve(struct fe *problem, matrix *guess)
{
    matrix *dx; /* How much to update the guess by */
    matrix *newguess;
    //int rows = (2*problem->mesh->nelemx+1)*(2*problem->mesh->nelemy+1);
    int rows = (problem->mesh->nelemx+1)*(problem->mesh->nelemy+1);
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
