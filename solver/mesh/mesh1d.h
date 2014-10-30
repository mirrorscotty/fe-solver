#ifndef MESH1D_H
#define MESH1D_H

#include "matrix.h"
#include "basis.h"

typedef struct {
    vector *points;
    vector *map;
} Elem1D;

typedef struct mesh1dstruct {
    int nelem; /**< Number of elements in the mesh */

    int nnodes; /**< Number of nodes */

    double x1; /**< Left-most node coordinate */
    double x2; /**< Right-most node coordinate */

    int t; /**< Earliest time index the mesh is valid for */
    struct mesh1dstruct *next; /**< Previous meshes for the domain */
    struct mesh1dstruct *orig; /**< Originally define mesh for the problem */

    Elem1D **elem; /**< Set of all of the elements in the mesh */
    vector *nodes; /**< Vector of node coordinates */
} Mesh1D;

#include "finite-element1d.h"

Mesh1D* GenerateUniformMesh1D(basis*, double, double, int);
void DestroyMesh1D(Mesh1D*);
Elem1D* CreateElem1D(basis*);
void DestroyElem1D(Elem1D*);
void meshprnt1d(Mesh1D*);

Mesh1D* Remesh1D(Mesh1D*, vector*);

#endif

