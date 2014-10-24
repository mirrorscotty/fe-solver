#ifndef SOLVE_H
#define SOLVE_H

#include "finite-element.h"
#include "finite-element1d.h"

/* 2D solvers */
matrix* LinSolve(struct fe*);
matrix* NLinSolve(struct fe*, matrix*);

/* 1D, Linear Solvers */
matrix* LinSolve1D(struct fe1d*);
matrix* LinSolve1DTrans(struct fe1d*);
matrix* LinSolve1DTransImp(struct fe1d*);

/* 1D, Nonlinear Solvers */
matrix* NLinSolve1D(struct fe1d*, matrix*);
matrix* NLinSolve1DTransImp(struct fe1d*, matrix *guess);

/* 1D Prediction Algorithms */
matrix* PredictSolnO0(struct fe1d*);
matrix* PredictSolnO1(struct fe1d*);
matrix* PredictSolnO2(struct fe1d*);

#endif

