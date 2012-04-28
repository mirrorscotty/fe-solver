#ifndef MTXSOLVER_H
#define MTXSOLVER_H

#include "finite-element.h"

void ForwardSubstitution(matrix*);
void ReverseElimination(matrix*);
matrix* SolveMatrixEquation(matrix*, matrix*);
matrix* LinSolve(struct fe*);
matrix* NLinSolve(struct fe*, matrix*);

#endif
