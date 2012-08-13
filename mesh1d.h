#ifndef MESH1D_H
#define MESH1D_H

#include "matrix.h"
#include "basis.h"

typedef struct {
    vector *points;
    vector *map;
} Elem1D;

typedef struct {
    int nelem;

    int nnodes;

    double x1;
    double x2;

    Elem1D **elem;
    vector *nodes;
} Mesh1D;

Mesh1D* GenerateUniformMesh1D(basis*, double, double, int);
void DestroyMesh1D(Mesh1D*);
Elem1D* CreateElem1D(basis*);
void DestroyElem1D(Elem1D*);
void meshprnt1d(Mesh1D*);

#endif

