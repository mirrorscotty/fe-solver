/**
 * @mesh1d.h
 * Defines the mesh and element structures for a one-dimensional finite element
 * problem
 */

#ifndef MESH1D_H
#define MESH1D_H

#include "matrix.h"
#include "basis.h"

/**
 * Defines a single element in the mesh
 */
typedef struct elem1dstruct {
    vector *points; /**< Vector of the x-values of the nodes in the element. */
    vector *map; /**< Vector of the mesh node numbers for each node in the
                  *   element. These show where each point in the element is in
                  *   the mesh. */
    struct elem1dstruct *prev; /**< If this mesh has been deformed, this will
                                *   point to the element at the previous time
                                *   step. */
} Elem1D;

/**
 * One-dimensional finite element mesh structure.
 */
typedef struct mesh1dstruct {
    int nelem; /**< Number of elements in the mesh */

    int nnodes; /**< Number of nodes */

    double x1; /**< Left-most node coordinate */
    double x2; /**< Right-most node coordinate */

    int t; /**< Earliest time index the mesh is valid for */
    struct mesh1dstruct *prev; /**< Previous meshes for the domain */
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
int MeshIsOrig(Mesh1D*);

#endif

