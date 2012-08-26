#ifndef SOLUTION_H
#define SOLUTION_H

#include "matrix.h"
#include "mesh1d.h"

struct fe1d;

typedef struct {
    int t;
    double dt;
    matrix *val; /* The values of the variable at the indicated time. */
    matrix *dval; /* The values of the first time derivative. */
} solution;

solution* CreateSolution(int, double, matrix*);
void DestroySolution(solution*);

double EvalSoln1D(struct fe1d*, int, Elem1D*, solution*, double);

#endif

