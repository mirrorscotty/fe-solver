#ifndef MESH_H
#define MESH_H

#include "matrix.h"

typedef struct {
    double dx;
    double dy;
} Elem2D;

typedef struct {
    int nelemx;
    int nelemy;
    
    double x1;
    double x2;
    double y1;
    double y2;

    Elem2D *elem;
    //double *dx, *dy;
} Mesh2D;

Mesh2D* GenerateUniformMesh2D(double, double, double, double, int, int);
void DestroyMesh2D(Mesh2D*);
vector* GetNodeCoordinates(Mesh2D*, int);

#endif
