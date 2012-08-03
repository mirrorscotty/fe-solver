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
    double displace;
    double tol; /* Tolerance for the nonlinear solver */
    
    matrix* (*makej)(struct fe*, Elem2D*, matrix*);
    matrix* (*makef)(struct fe*, Elem2D*, matrix*);
    
    void (*applybcs)(struct fe*);
    
    /* Problem-specific parameters */
    double P;
    double a;
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
    
#endif
