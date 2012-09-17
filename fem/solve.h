#ifndef SOLVE_H
#define SOLVE_H

#include "finite-element.h"
#include "finite-element1d.h"

matrix* LinSolve(struct fe*);
matrix* LinSolve1D(struct fe1d*);
matrix* LinSolve1DTrans(struct fe1d*);
matrix* LinSolve1DTransImp(struct fe1d*);
matrix* NLinSolve(struct fe*, matrix*);

#endif
