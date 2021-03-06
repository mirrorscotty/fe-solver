#ifndef MESH_H
#define MESH_H

#include "matrix.h"
#include "solver/integration/basis.h"

typedef struct {
    //double dx; /* TODO: Remove. Should no longer be used anywhere */
    //double dy; /* TODO: Remove. Should no longer be used anywhere */
    vector** points;
    vector* map;
} Elem2D;

/* For 1D meshes, just use the x-coordinate. The variable nelemy should be 1 */
typedef struct {
    int nelemx;
    int nelemy;

    int nnodes; /* Total number of nodes in the mesh */
    /* Todo: Use this variable ^^ */

    double x1;
    double x2;
    double y1;
    double y2;

    Elem2D **elem;
    vector **nodes;
    //double *dx, *dy;
} Mesh2D;

Mesh2D* GenerateUniformMesh2D(double, double, double, double, int, int);
Mesh2D* GenerateRadialMesh(basis*, vector*, vector*, double, int, int);
vector* MakeR(basis*, double, double, double, int);
Mesh2D* MakeSpheroidMesh(basis*, double, double, int, int);

void DestroyMesh2D(Mesh2D*);
void MeshPrint(Mesh2D*);
void MeshPrintNodes(Mesh2D*);
vector* GetNodeCoordinates(Mesh2D*, int);
Elem2D* CreateElem2D(basis*);
void DestroyElem2D(Elem2D*);

#endif
