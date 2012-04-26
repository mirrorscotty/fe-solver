#ifndef MESH_H
#define MESH_H

#include "matrix.h"

typedef struct {
    double dx; /* TODO: Remove. Should no longer be used anywhere */
    double dy; /* TODO: Remove. Should no longer be used anywhere */
    vector** points;
    vector* map;
} Elem2D;

typedef struct {
    int nelemx;
    int nelemy;
    
    double x1;
    double x2;
    double y1;
    double y2;

    Elem2D **elem;
    //double *dx, *dy;
} Mesh2D;

Mesh2D* GenerateUniformMesh2D(double, double, double, double, int, int);
Mesh2D* GenerateRadialMesh(vector*, vector*, double, int, int);
vector* MakeR(double, double, double, int);
Mesh2D* MakeSpheroidMesh(double, double, int, int);

void DestroyMesh2D(Mesh2D*);
vector* GetNodeCoordinates(Mesh2D*, int);
Elem2D* CreateElem2D();
void DestroyElem2D(Elem2D*);

#endif
