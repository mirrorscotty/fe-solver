#ifndef MTXSOLVER_H
#define MTXSOLVER_H

#include "finite-element.h"
#include "finite-element1d.h"

void ForwardSubstitution(matrix*);
void ReverseElimination(matrix*);
matrix* SolveMatrixEquation(matrix*, matrix*);
matrix* LinSolve(struct fe*);
matrix* LinSolve1D(struct fe1d*);
matrix* LinSolve1DTrans(struct fe1d*);
matrix* LinSolve1DTransImp(struct fe1d*);
matrix* NLinSolve(struct fe*, matrix*);

#endif
