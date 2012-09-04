#ifndef FINITE_ELEMENT_H
#define FINITE_ELEMENT_H

#include "matrix.h"
#include "basis.h"
#include "mesh2d.h"

struct fe;

struct fe {
    basis *b;
    Mesh2D *mesh;
    
    matrix *J;
    //matrix *Jconst; /* Jacobian with constraints. If there are no constraints,
    //                   then this points to the same thing as J. */
    matrix *F;
    matrix *R;
    
    int nvars; /* Number of dependant variables. */
    int nrows; /* Number of rows (nodes) in the final solution. */
    int nconstr; /* Number of additional constraints (arc length, etc.) */

    double tol; /* Tolerance for the nonlinear solver */
    
    /* Functions to make the element mass/stiffness matrix and load vector */
    matrix* (*makej)(struct fe*, Elem2D*, matrix*);
    matrix* (*makef)(struct fe*, Elem2D*, matrix*);
    
    /* Function that applies boundary conditions to the global matricies */
    void (*applybcs)(struct fe*);
    
    /* Problem-specific parameters */
    double P; /* Pressure applied to the top of the column */
    double a; /* Ratio of pressure on the top surface to the side surface */
    //double displace; /* Displacement of the wall for contact problems */
    
    /* Arc length to use when using the nonlinear solver */
    double ds;
    
    /* Pointers to the current and previous guess (for nonlinear solver) */
    /* These aren't currently used. */
    matrix *guess;
    matrix *prevguess;
};

struct fe* CreateFE(basis*,
                    Mesh2D*,
                    matrix* (*)(struct fe*, Elem2D*, matrix*),
                    matrix* (*)(struct fe*, Elem2D*, matrix*),
                    void (*)(struct fe*));

void DestroyFE(struct fe*);

matrix* AssembleJ(struct fe*, matrix*);
matrix* AssembleF(struct fe*, matrix*);
matrix* CalcResidual(struct fe*, matrix*);
matrix* GetLocalGuess(struct fe*, matrix*, int);
void ApplyEssentialBC(struct fe*, int, int (*)(struct fe*, int), double (*)(struct fe*, int));
void ApplyNaturalBC(struct fe*, int, int (*)(struct fe*, int), double (*)(struct fe*, int));
#endif

