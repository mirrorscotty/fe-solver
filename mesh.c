#include <stdlib.h>
#include <stdio.h>
#include "matrix.h"
#include "mesh.h"

Mesh2D* GenerateUniformMesh2D(double x1, double x2,
                              double y1, double y2,
                              int nx, int ny)
{
    int i, j, z;
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

    mesh->elem = (Elem2D*) calloc(nx*ny, sizeof(Elem2D*));

    /* TODO: automatically change the direction in which the elements are
     * numbered for maximum efficiency. */
    for(i=0; i<ny; i++) {
        for(j=0; j<nx; j++) {
            z = i*ny+j; // Abreviate the the index.
            mesh->elem[z] = CreateElem2D();

            setvalV(mesh->elem[z]->points[0], 0, x1+j*nx);
            setvalV(mesh->elem[z]->points[0], 1, y1+i*ny);

            setvalV(mesh->elem[z]->points[1], 0, x1+j*nx);
            setvalV(mesh->elem[z]->points[1], 1, y1+(i+1)*ny);
            
            setvalV(mesh->elem[z]->points[2], 0, x1+(j+1)*nx);
            setvalV(mesh->elem[z]->points[2], 1, y1+i*ny);

            setvalV(mesh->elem[z]->points[3], 0, x1+(j+1)*nx);
            setvalV(mesh->elem[z]->points[3], 1, y1+(i+1)*ny);
        }
    }

    return mesh;
}

vector* GetNodeCoordinates(Mesh2D *mesh, int node)
{
    vector *v;
    v = CreateVector(2);

    int x, y;
    double dx, dy;

    x = node / (mesh->nelemy+1);
    y = node - x*(mesh->nelemy+1);

    // Cheat because the mesh is always uniform.
    //dx = mesh->elem[0].dx;
    //dy = mesh->elem[0].dy;

//    printf("Node: %g, %g\n", dx*x, dy*y);

    setvalV(v, 0, x*dx);
    setvalV(v, 1, y*dy);

//    PrintVector(v);

    return v;
}

Elem2D* CreateElem2D()
{
    Elem2D *elem;
    int i;
    elem = (Elem2D*) calloc(1, sizeof(Elem2D));

    elem->points = (vector**) calloc(4, sizeof(vector*));
    for(i=0; i<4; i++)
        elem->points[i] = CreateVector(2);

    return elem;
}

void DestroyElem2D(Elem2D *elem)
{
    int i;
    // 4 nodes per element.
    for(i=0; i<4; i++)
        DestroyVector(elem->points[i]);
    free(elem->points);
}
        
void DestroyMesh2D(Mesh2D *mesh)
{
    free(mesh->elem);
    free(mesh);
}

