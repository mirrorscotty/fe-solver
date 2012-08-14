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
    s->values = values;

    return s;
}

void DestroySolution(solution *s)
{
    DestroyMatrix(s->values);
    free(s);
}

/* This function needs to be modified to support hermite cubic basis functions.
 */
double IEvalSoln1D(struct fe1d *p, Elem1D *elem, solution *s, double xi)
{   
    int i;
    double result = 0;
    int n = p->b->n;

    for(i=0; i<n; i++) {
        result += p->b->phi[i](xi) * val(s->values, valV(elem->map, i), 0);
    }
    
    return result;
} 

