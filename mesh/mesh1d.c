#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "matrix.h"
#include "mesh1d.h"

void meshprnt1d(Mesh1D* mesh)
{
    int i;
    printf("Elements: %d\n", mesh->nelem);
    printf("---------------------\n");
    for(i=0; i<mesh->nelem; i++) {
        printf("Element #%d\n---------------------\n", i);
        PrintVector(mesh->elem[i]->map);
        PrintVector(mesh->elem[i]->points);
        puts("");
    }
    return;
}

Mesh1D* GenerateUniformMesh1D(basis *b, double x1, double x2, int nx)
{
    int i, j;
    double dx = (x2-x1)/nx;

    Mesh1D *mesh;
    mesh = (Mesh1D*) calloc(1, sizeof(Mesh1D));

    mesh->x1 = x1;
    mesh->x2 = x2;

    mesh->nelem = nx;

    mesh->elem = (Elem1D**) calloc(nx, sizeof(Elem1D));
    mesh->nodes = CreateVector(nx+1);

    for(i=0; i<nx; i++) {
        mesh->elem[i] = CreateElem1D(b);
        setvalV(mesh->elem[i]->points, 0, x1+dx*i);
        setvalV(mesh->elem[i]->points, 1, x1+dx*(i+1));

        setvalV(mesh->elem[i]->map, 0, i);
        setvalV(mesh->elem[i]->map, 1, i+1);

        for(j=0; j<2; j++) {
            setvalV(mesh->nodes, (int) valV(mesh->elem[i]->map, j),
                    valV(mesh->elem[i]->points, j));
        }
    }

    return mesh;
}

Elem1D* CreateElem1D(basis *b)
{
    Elem1D *elem;
    int nnodes = b->n;
    elem = (Elem1D*) calloc(1, sizeof(Elem1D));
    elem->points = CreateVector(nnodes);
    elem->map = CreateVector(nnodes);

    return elem;
}

void DestroyElem1D(Elem1D *elem)
{
    DestroyVector(elem->points);
    DestroyVector(elem->map);
    free(elem);
}

void DestroyMesh1D(Mesh1D *mesh)
{
    int i;
    for(i=0; i<mesh->nelem; i++)
        DestroyElem1D(mesh->elem[i]);
    DestroyVector(mesh->nodes);
    free(mesh->elem);
    free(mesh);
}

