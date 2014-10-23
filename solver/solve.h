#ifndef SOLVE_H
#define SOLVE_H

#include "finite-element.h"
#include "finite-element1d.h"

matrix* LinSolve(struct fe*);
matrix* LinSolve1D(struct fe1d*);
matrix* LinSolve1DTrans(struct fe1d*);
matrix* LinSolve1DTransImp(struct fe1d*);
matrix* NLinSolve1DTransImp(struct fe1d*, matrix *guess);

matrix* NLinSolve(struct fe*, matrix*);
matrix* NLinSolve1D(struct fe1d*, matrix*);

matrix* CalcTimeDerivative(struct fe1d*, matrix*);

matrix* PredictSolnO0(struct fe1d*);
matrix* PredictSolnO1(struct fe1d*);
matrix* PredictSolnO2(struct fe1d*);

#endif

