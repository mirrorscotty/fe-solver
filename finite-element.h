#ifndef FINITE_ELEMENT_H
#define FINITE_ELEMENT_H

#include "matrix.h"
#include "basis.h"
#include "mesh.h"

struct fe;

struct fe {
    basis *b;
    Mesh2D *mesh;
    
    matrix *J;
    matrix *F;
    matrix *R;
    
    int nvars; /* Number of variables being solved for. */
    double tol; /* Tolerance for the nonlinear solver */
    
    matrix* (*makej)(struct fe*, matrix*);
    matrix* (*makef)(struct fe*, matrix*);
    
    void (*applybcs)(struct fe*);
    
    /* Problem-specific parameters */
    double P;
    double a;
};

struct fe* CreateFE(basis*,
                    Mesh2D*,
                    matrix* (*)(struct fe*, matrix*),
                    matrix* (*)(struct fe*, matrix*),
                    void (*)(struct fe*));

void DestroyFE(struct fe*);

matrix* AssembleJ(struct fe*, matrix*);
matrix* AssembleF(struct fe*, matrix*);
matrix* CalcResidual(struct fe*, matrix*);
matrix* GetLocalGuess(struct fe*, matrix*, int);
    
#endif