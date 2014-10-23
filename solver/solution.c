#include <stdlib.h>

#include "solution.h"
#include "matrix.h"

solution* CreateSolution(int tindex, double deltat, matrix *values)
{
    solution *s;
    s = (solution*) calloc(1, sizeof(solution));

    s->t = tindex;
    s->dt = deltat;
    s->val = values;
    s->dval = NULL;

    return s;
}

void DestroySolution(solution *s)
{
    if(s) {
        if(s->val)
            DestroyMatrix(s->val);
        if(s->dval)
            DestroyMatrix(s->dval);
        free(s);
    }
}

void DeleteTimeDeriv(solution *s) {
    //if(s->dval) {
        DestroyMatrix(s->dval);
        s->dval = NULL;
    //}
}

