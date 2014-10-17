#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "matrix.h"
#include "mesh1d.h"

/**
 * @brief Print out all of the points in a 1D mesh
 * @param mesh The mesh to print
 */
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

/**
 * Create a one-dimensional mesh with uniformly spaced nodes.
 * @param b Set of basis functions to use
 * @param x1 Coordinate of the left-most node
 * @param x2 Coordinate of the right-most node
 * @param nx Number of nodes to put into the mesh
 * @returns 1D mesh structure
 */
Mesh1D* GenerateUniformMesh1D(basis *b, double x1, double x2, int nx)
{
    int i, j, k;
    double dx = (x2-x1)/nx;
    int n = b->n; /* Number of nodes per element */

    Mesh1D *mesh;
    mesh = (Mesh1D*) calloc(1, sizeof(Mesh1D));

    mesh->x1 = x1;
    mesh->x2 = x2;

    mesh->nelem = nx;

    mesh->elem = (Elem1D**) calloc(nx, sizeof(Elem1D));
    mesh->nodes = CreateVector((n-1)*nx+1);

    for(i=0; i<nx; i++) {
        mesh->elem[i] = CreateElem1D(b);
        for(k=0; k<n; k++) {
            setvalV(mesh->elem[i]->points, k, x1+dx*(i*(n-1)+k)/(n-1));

            setvalV(mesh->elem[i]->map, k, i*(n-1)+k);
        }

        for(j=0; j<n; j++) {
            setvalV(mesh->nodes, (int) valV(mesh->elem[i]->map, j),
                    valV(mesh->elem[i]->points, j));
        }
    }

    return mesh;
}

/**
 * Copy a mesh and all of the elements in it. This is useful for stuff like
 * remeshing or allowing the mesh to move/shrink.
 * @param orig Pointer to the original mesh
 * @returns A pointer to a copy of the original
 */
Mesh1D* CopyMesh1D(Mesh1D* orig)
{
    int i;

    Mesh1D *copy;
    copy = (Mesh1D*) calloc(1, sizeof(Mesh1D));
    copy->elem = (Elem1D**) calloc(orig->nelem, sizeof(Elem1D));
    copy->nodes = CopyVector(orig->nodes);

    copy->x1 = orig->x1;
    copy->x2 = orig->x2;
    copy->nelem = orig->nelem;
    /* Don't copy nnodes because nothing uses it */
    copy->t = orig->t;
    copy->next = NULL; /* Set the next node to NULL */

    for(i=0; i<orig->nelem; i++) {
        copy->elem[i] = (Elem1D*) calloc(1, sizeof(Elem1D));
        copy->elem[i]->points = CopyVector(orig->elem[i]->points);
        copy->elem[i]->map = CopyVector(orig->elem[i]->map);
    }

    return copy;
}

/**
 * Make a new mesh where the nodes have moved from their locations in the 
 * original mesh to the locations supplied.
 * @param orig Mesh to get the element information, etc. from
 * @param nodes Vector of new node locations
 * @returns A mesh with the same number of points as the original with new
 *      coordinates.
 */
Mesh1D* Remesh1D(Mesh1D* orig, vector* nodes)
{
    Mesh1D *new;
    int i, j, k;

    new = CopyMesh1D(orig);
    DestroyVector(new->nodes);
    new->nodes = nodes;

    for(i=0; i<nelem; i++) {
        for(j=0; j<len(new->elem[i]->points); j++) {
            setvalV(new->elem[i]->points,
                    j,
                    valV(new->nodes, valV(new->elem[i]->map, j)));
        }
    }

    return new;
}

Mesh1D* MoveMeshF(struct fe1d *p, Mesh1D *orig, double t
                  double (*F)(struct fe1d *p, double x, double t))
{
    int i;
    double dx;
    vector *new;

    new = CreateVector(len(orig->nodes));

    setvalV(new, 0, valV(orig->nodes, 0)); /* The first node doesn't change */

    for(i=1; i<len(orig->nodes); i++) {
        dx = valV(orig->nodes, i) - valV(orig->nodes, i-1);
        x = valV(orig->nodes, i);
        F = F(p, x, t);
        setvalV(new, i, F*dx+valV(new, i-1));
    }

    return Remesh(orig, new);
}

/**
 * Create a single element to insert into a mesh. These get passed to the
 * integration functions so that they know the points they should be using to
 * integrate over. This structure holds all of the points for a single element
 * in a vector.
 * @param b Set of basis functions to use when making the node.
 * @returns 1D element structure
 */
Elem1D* CreateElem1D(basis *b)
{
    Elem1D *elem;
    int nnodes = b->n;
    elem = (Elem1D*) calloc(1, sizeof(Elem1D));
    elem->points = CreateVector(nnodes);
    elem->map = CreateVector(nnodes);

    return elem;
}

/**
 * Delete a single mesh element
 * @param elem Pointer to the element to delete
 */
void DestroyElem1D(Elem1D *elem)
{
    DestroyVector(elem->points);
    DestroyVector(elem->map);
    free(elem);
}

/**
 * Delete a mesh and all associated elements.
 * @param mesh Pointer to the mesh to delete
 */
void DestroyMesh1D(Mesh1D *mesh)
{
    int i;
    for(i=0; i<mesh->nelem; i++)
        DestroyElem1D(mesh->elem[i]);
    DestroyVector(mesh->nodes);
    free(mesh->elem);
    free(mesh);
}

