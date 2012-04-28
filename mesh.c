#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "matrix.h"
#include "mesh.h"

/* Print out a mesh in a somewhat readable format */
void MeshPrint(Mesh2D *mesh)
{
    int i, j;
    printf("Elements: %d (%d, %d)\n", mesh->nelemx*mesh->nelemy,
            mesh->nelemx, mesh->nelemy);
    printf("---------------------\n");
    for(i=0; i<mesh->nelemx*mesh->nelemy; i++) {
        printf("Element #%d\n---------------------\n", i);
        for(j=0; j<len(mesh->elem[i]->map); j++) {
            PrintVector(mesh->elem[i]->points[j]);
        }
        puts("");
    }
    return;
}

/* Generate a uniform mesh for a spheroid. r0 is the vector of r values that
 * specifies the size and shape of the hole in the center of the mesh, and r1 is
 * the set of r values defining the outside edge of the mesh. The number of
 * values in r0 and r1 should be equal to Ntheta (the number of nodes in the
 * theta direction). The variable "thetamax" describes the portion of the
 * spheroid which is meshed, and Nr and Ntheta are the number of elements in the
 * r and theta directions, respectively.
 *
 * This whole function needs to be de-jankified. The whole multiplying i and j
 * by two needs to go and the sqrt(b->n) should disappear as well.
 */
Mesh2D* GenerateRadialMesh(basis *b, vector *r0, vector *r1, double thetamax, 
                           int Nr, int Ntheta)
{
    int i, j, k, z;
    int c, d;
    double ri, ri1, ri1j1, rij1, thetai, thetai1;
    /* Just allocate enough for a BiQuad mesh since it's easy */
    double r[3][3], theta[3]; 
    int nnodes = b->n; /* Number of nodes per element */

    Elem2D *e;
    Mesh2D *mesh;
    mesh = (Mesh2D*) calloc(1, sizeof(Mesh2D));

    /* This isn't exactly correct, but it's close enough. */
    mesh->nelemx = Ntheta;
    mesh->nelemy = Nr;

    mesh->elem = (Elem2D**) calloc(Nr*Ntheta, sizeof(Elem2D*));
    mesh->nodes = (vector**) calloc((Nr+1)*(Ntheta+1), sizeof(vector*));

    for(i=0; i<Nr; i++) {
        for(j=0; j<Ntheta; j++) {
            z = i*Nr+j; /* Element number */
            printf("%d\n", z);
            mesh->elem[z] = CreateElem2D(b);
            e = mesh->elem[z];

            for(c=0; c<sqrt(nnodes); c++) {
                for(d=0; d<sqrt(nnodes); d++) {
                    r[c][d] = (valV(r1, 2*j+d) - valV(r0, 2*j+d))/Nr/(sqrt(b->n)-b->overlap) * (2*i+c) + valV(r0, 2*j+d);
                    theta[d] = thetamax/Ntheta/(sqrt(b->n)-b->overlap) * (2*j+d);
                }
            }

            for(c=0; c<sqrt(nnodes); c++) {
                for(d=0; d<sqrt(nnodes); d++) {
                    setvalV(e->points[(int) (sqrt(nnodes)*d+c)], 0, r[c][d]*cos(theta[d]));
                    setvalV(e->points[(int) (sqrt(nnodes)*d+c)], 1, r[c][d]*sin(theta[d]));
                }
            }

            /* Determine the global node numbers for each node in the
             * element. Used for matrix assembly. */
            /*
            setvalV(e->map, 0, (double) i*(Nr+1)+j);
            setvalV(e->map, 1, (double) (i+1)*(Nr+1)+j);
            setvalV(e->map, 2, (double) i*(Nr+1)+j+1);
            setvalV(e->map, 3, (double) (i+1)*(Nr+1)+j+1);
            */

            /*
            for(k=0;k<4;k++) {
                mesh->nodes[(int) valV(e->map, k)] = e->points[k];
            }
            */
        }
    }

    //MeshPrint(mesh);
    return mesh;
}

/* Make the R vectors for so that a radial mesh can be constructed. */
vector* MakeR(basis *bas, double a, double b, double thetamax, int Ntheta)
{
    int i;
    double theta;
    vector *result;

    result = CreateVector(Ntheta+bas->n-bas->overlap);

    for(i=0; i<=Ntheta+bas->n-bas->overlap; i++) {
        theta = thetamax/Ntheta*i;
        setvalV(result, i, 1/sqrt(pow(cos(theta), 2)/(a*a) + pow(sin(theta), 2)/(b*b)));
    }


    return result;
}

Mesh2D* MakeSpheroidMesh(basis *bas, double e, double r, int Nr, int Nt)
{
    Mesh2D *mesh;
    vector *r0, *r1;
    double a, b;
    a = 1; /* Set the major axis of the spheroid to 1 */
    b = sqrt((1-e*e)*a*a); /* Calculate the minor axis from a and the 
                            * eccentricity */
    r0 = MakeR(bas, b, a, M_PI_2, Nt);
    r1 = MakeR(bas, r, r, M_PI_2, Nt);

    mesh = GenerateRadialMesh(bas, r0, r1, M_PI_2, Nr, Nt);

    DestroyVector(r0);
    DestroyVector(r1);

    return mesh;
}

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

    mesh->elem = (Elem2D**) calloc(nx*ny, sizeof(Elem2D*));

    /* TODO: automatically change the direction in which the elements are
     * numbered for maximum efficiency. */
    for(i=0; i<ny; i++) {
        for(j=0; j<nx; j++) {
            z = i*ny+j; // Abreviate the the index.
            mesh->elem[z] = CreateElem2D(NULL); /* DOESN"T WORK! */

            setvalV(mesh->elem[z]->points[0], 0, x1+j*DeltaX);
            setvalV(mesh->elem[z]->points[0], 1, y1+i*DeltaY);

            setvalV(mesh->elem[z]->points[1], 0, x1+j*DeltaX);
            setvalV(mesh->elem[z]->points[1], 1, y1+(i+1)*DeltaY);
            
            setvalV(mesh->elem[z]->points[2], 0, x1+(j+1)*DeltaX);
            setvalV(mesh->elem[z]->points[2], 1, y1+i*DeltaY);

            setvalV(mesh->elem[z]->points[3], 0, x1+(j+1)*DeltaX);
            setvalV(mesh->elem[z]->points[3], 1, y1+(i+1)*DeltaY);

            /* TODO: Add in node numbers */

        }
    }

    return mesh;
}

/* TODO: Remove this function. Or, at the very least, rewrite it so that it
 * makes sense. */
vector* GetNodeCoordinates(Mesh2D *mesh, int node)
{
    vector *v;
    v = CreateVector(2);

    int x, y;
    double dx = 0, dy = 0;

    x = node / (mesh->nelemy+1);
    y = node - x*(mesh->nelemy+1);

    // Cheat because the mesh is always uniform.
    //dx = mesh->elem[0].dx;
    //dy = mesh->elem[0].dy;

    setvalV(v, 0, x*dx);
    setvalV(v, 1, y*dy);

//    PrintVector(v);

    return v;
}


Elem2D* CreateElem2D(basis *b)
{
    Elem2D *elem;
    int i;
    int nnodes = b->n;
    elem = (Elem2D*) calloc(1, sizeof(Elem2D));

    elem->points = (vector**) calloc(nnodes, sizeof(vector*));
    for(i=0; i<nnodes; i++)
        elem->points[i] = CreateVector(2);

    elem->map = CreateVector(nnodes);

    return elem;
}

void DestroyElem2D(Elem2D *elem)
{
    int i;
    // The length of the map is equal to the number of nodes
    for(i=0; i<len(elem->map); i++)
        DestroyVector(elem->points[i]);
    free(elem->points);
    DestroyVector(elem->map);

    free(elem);
}
        
void DestroyMesh2D(Mesh2D *mesh)
{
    int i;
    for(i=0; i<mesh->nelemx*mesh->nelemy; i++)
        free(mesh->elem[i]);
    free(mesh->nodes);
    free(mesh->elem);
    free(mesh);
}

