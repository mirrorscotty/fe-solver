#include <stdlib.h>

#include "solution.h"
#include "matrix.h"
#include "finite-element1d.h"

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

/* This function needs to be modified to support hermite cubic basis functions.
 */
double EvalSoln1D(struct fe1d *p, int var, Elem1D *elem, solution *s, double xi)
{   
    int i;
    double result = 0;
    int n = p->b->n;
    int nvars = p->nvars;

    for(i=0; i<n; i++) {
        result += p->b->phi[i](xi)
                  * val(s->val, valV(elem->map, i)*nvars, 0);
    }
    
    return result;
} 

void DeleteTimeDeriv(solution *s) {
    //if(s->dval) {
        DestroyMatrix(s->dval);
        s->dval = NULL;
    //}
}
