/**
 * @file solution.c
 * Set of functions for managing PDE solutions
 */

#include <stdlib.h>

#include "solution.h"
#include "matrix.h"

/**
 * Create a solution structure and store a set of dependent variable values
 * in it. The time derivative of these values is not specified.
 * @param tindex Time step number
 * @param deltat Time step size
 * @param values Column matrix of dependent variable values. If multiple PDEs
 *      are being solved, then the matrix will have N*v rows, where N is the
 *      number of nodes in the mesh, and v is the number of PDEs. In this case,
 *      all the dependent variable values at the first node will be first, then
 *      the values at the second node (in the same order), and so on.
 * @returns A solution structure with all the values stored in it.
 */
solution* CreateSolution(int tindex, double deltat, matrix *values)
{
    solution *s;
    /* Allocate memory */
    s = (solution*) calloc(1, sizeof(solution));

    /* Store the values */
    s->t = tindex;
    s->dt = deltat;
    s->val = values;
    s->dval = NULL;

    return s;
}

/**
 * Deallocate memory for a solution structure. This will also delete any
 * matricies stored in it.
 * @param s Solution to delete
 */
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

/**
 * This deletes the stored time derivative (mostly to save memory).
 * @param s Solution structure
 */
void DeleteTimeDeriv(solution *s)
{
    //if(s->dval) {
        DestroyMatrix(s->dval);
        s->dval = NULL;
    //}
}

