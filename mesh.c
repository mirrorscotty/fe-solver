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
        for(j=0; j<4; j++) {
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
 */
Mesh2D* GenerateRadialMesh(vector *r0, vector *r1, double thetamax, 
                           int Nr, int Ntheta)
{
    int i, j, z;
    double ri, ri1, ri1j1, rij1, thetai, thetai1;

    Elem2D *e;
    Mesh2D *mesh;
    mesh = (Mesh2D*) calloc(1, sizeof(Mesh2D));

    /* This isn't exactly correct, but it's close enough. */
    mesh->nelemx = Ntheta;
    mesh->nelemy = Nr;

    mesh->elem = (Elem2D**) calloc(Nr*Ntheta, sizeof(Elem2D*));

    for(i=0; i<Nr; i++) {
        for(j=0; j<Ntheta; j++) {
            z = i*Nr+j; /* Element number */
            mesh->elem[z] = CreateElem2D();
            e = mesh->elem[z];

            ri = (valV(r1, j) - valV(r0, j))/Nr * i + valV(r0, j);
            ri1 = (valV(r1, j) - valV(r0, j))/Nr * (i+1) + valV(r0, j);
            rij1 = (valV(r1, j+1) - valV(r0, j+1))/Nr * (i) + valV(r0, j+1);
            ri1j1 = (valV(r1, j+1) - valV(r0, j+1))/Nr * (i+1) + valV(r0, j+1);

            thetai = thetamax/Ntheta * j;
            thetai1 = thetamax/Ntheta * (j+1);

            setvalV(e->points[0], 0, ri*cos(thetai));
            setvalV(e->points[0], 1, ri*sin(thetai));

            setvalV(e->points[1], 0, ri1*cos(thetai));
            setvalV(e->points[1], 1, ri1*sin(thetai));

            setvalV(e->points[2], 0, rij1*cos(thetai1));
            setvalV(e->points[2], 1, rij1*sin(thetai1));

            setvalV(e->points[3], 0, ri1j1*cos(thetai1));
            setvalV(e->points[3], 1, ri1j1*sin(thetai1));
        }
    }

    return mesh;
}

/* Make the R vectors for so that a radial mesh can be constructed. */
vector* MakeR(double a, double b, double thetamax, int Ntheta)
{
    int i;
    double theta;
    vector *result;

    result = CreateVector(Ntheta+1);

    for(i=0; i<=Ntheta; i++) {
        theta = thetamax/Ntheta*i;
        setvalV(result, i, 1/sqrt(pow(cos(theta), 2)/(a*a) + pow(sin(theta), 2)/(b*b)));
    }


    return result;
}

Mesh2D* MakeSpheroidMesh(double e, double r, int Nr, int Nt)
{
    Mesh2D *mesh;
    vector *r0, *r1;
    double a, b;
    a = 1; /* Set the major axis of the spheroid to 1 */
    b = sqrt((1-e*e)*a*a); /* Calculate the minor axis from a and the 
                            * eccentricity */
    r0 = MakeR(b, a, M_PI_2, Nt);
    r1 = MakeR(r, r, M_PI_2, Nt);

    mesh = GenerateRadialMesh(r0, r1, M_PI_2, Nr, Nt);

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
            mesh->elem[z] = CreateElem2D();

            setvalV(mesh->elem[z]->points[0], 0, x1+j*DeltaX);
            setvalV(mesh->elem[z]->points[0], 1, y1+i*DeltaY);

            setvalV(mesh->elem[z]->points[1], 0, x1+j*DeltaX);
            setvalV(mesh->elem[z]->points[1], 1, y1+(i+1)*DeltaY);
            
            setvalV(mesh->elem[z]->points[2], 0, x1+(j+1)*DeltaX);
            setvalV(mesh->elem[z]->points[2], 1, y1+i*DeltaY);

            setvalV(mesh->elem[z]->points[3], 0, x1+(j+1)*DeltaX);
            setvalV(mesh->elem[z]->points[3], 1, y1+(i+1)*DeltaY);
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

