#ifndef FINITE_ELEMENT1D_H
#define FINITE_ELEMENT1D_h

#include "matrix.h"
#include "basis.h"
#include "mesh1d.h"

struct fe1d;

struct fe1d {
    basis *b;
    Mesh1D *mesh;

    matrix *J;
    matrix *R;
    matrix *F;

    int nvars; /* Number of dependant variables. */
    int nrows; /* Number of rows (nodes) in the final solution. */

    double tol; /* Tolerance for the nonlinear solver */

    /* Functions to make the element mass/stiffness matrix and the load vector */
    matrix* (*makej)(struct fe1d*, Elem1D*, matrix*);
    matrix* (*makek)(struct fe1d*, Elem1D*, matrix*);

    /* Function to apply all the relevant boundary conditions for the problem. */
    void (*applybcs)(struct fe1d*);
}

struct fe1d* CreateFE1D(basis*, Mesh1D*,
                        matrix* (*)(struct fe1d*, Elem1D*, matrix*),
                        matrix* (*)(struct fe1d*, Elem1D*, matrix*),
                        void (*)(struct fe1d*));
void DestroyFE1D(struct fe1d*);

matrix* AssembleJ1D(struct fe1d*, matrix*);
matrix* AssembleF1D(struct fe1d*, matrix*);

