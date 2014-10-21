#ifndef SOLUTION_H
#define SOLUTION_H

#include "matrix.h"

typedef struct {
    int t;
    double dt;
    matrix *val; /* The values of the variable at the indicated time. */
    matrix *dval; /* The values of the first time derivative. */
} solution;

solution* CreateSolution(int, double, matrix*);
void DestroySolution(solution*);
void DeleteTimeDeriv(solution*);


#endif

