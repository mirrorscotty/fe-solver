#include <stdlib.h>
#include <stdio.h>
#include "matrix.h"
#include "mesh.h"

Mesh2D* GenerateUniformMesh2D(double x1, double x2,
                              double y1, double y2,
                              int nx, int ny)
{
    int i;
    double DeltaX, DeltaY;
    DeltaX = (x2-x1)/nx;
    DeltaY = (y2-y1)/ny;

    Mesh2D *mesh;
    mesh = (Mesh2D*) calloc(1, sizeof(Mesh2D));

    mesh->x1 = x1;
    mesh->x2 = x2;
    mesh->y1 = y1;
    mesh->y2 = y2;

    mesh->nelemx = nx;
    mesh->nelemy = ny;

    mesh->elem = (Elem2D*) calloc(nx*ny, sizeof(Elem2D));
    
    for(i=0; i<nx*ny; i++) {
        mesh->elem[i].dx = DeltaX;
        mesh->elem[i].dy = DeltaY;
    }

    return mesh;
}

vector* GetNodeCoordinates(Mesh2D *mesh, int node)
{
    vector *v;
    v = CreateVector(2);

    int i, x, y;
    double dx, dy;

    x = node / (mesh->nelemy+1);
    y = node - x*(mesh->nelemy+1);

    // Cheat because the mesh is always uniform.
    dx = mesh->elem[0].dx;
    dy = mesh->elem[0].dy;

//    printf("Node: %g, %g\n", dx*x, dy*y);

    setvalV(v, 0, x*dx);
    setvalV(v, 1, y*dy);

//    PrintVector(v);

    return v;
}
        
void DestroyMesh2D(Mesh2D *mesh)
{
    free(mesh->elem);
    free(mesh);
}

