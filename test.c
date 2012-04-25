/*
 * Series of test functions for figuring stuff out. Pretty much nothing here is
 * used for actual analysis.
 */

#include <stdio.h>
#include "basis.h"
#include "mesh.h"
#include "isoparam.h"
#include "finite-element.h"

int main(int argc, char *argv[])
{
    Mesh2D *mesh;
    struct fe *p;
    int i = 5;

    //mesh = GenerateUniformMesh2D(0, 10, 0, 10, 10, 10);
    mesh = MakeSpheroidMesh(.8, 5, 5, 5);

    p = CreateFE(MakeLinBasis(2), mesh, NULL, NULL, NULL); 

    PrintVector(mesh->elem[i]->points[0]);
    PrintVector(mesh->elem[i]->points[1]);
    PrintVector(mesh->elem[i]->points[2]);
    PrintVector(mesh->elem[i]->points[3]);

    printf("dx/dxi = %g\n", IMapXXi(p, mesh->elem[i], .5, .5));
    printf("dx/deta = %g\n", IMapXEta(p, mesh->elem[i], .5, .5));
    printf("dy/dxi = %g\n", IMapYXi(p, mesh->elem[i], .5, .5));
    printf("dy/deta = %g\n", IMapYEta(p, mesh->elem[i], .5, .5));

    //MeshPrint(mesh);

    DestroyMesh2D(mesh);
    return 0;
}

/* Verify that the node locations are where they should be. Just prints out the
 * values of the four basis functions for a bilinear element at each of the four
 * nodes. All the values should be 1.
 */
void PrintNodeLocations(basis *b)
{
    printf("Phi0(%d, %d) = %g\n", 0, 0, EvalLin2D(b, 0, 0, 0));
    printf("Phi1(%d, %d) = %g\n", 0, 1, EvalLin2D(b, 1, 0, 1));
    printf("Phi2(%d, %d) = %g\n", 1, 0, EvalLin2D(b, 2, 1, 0));
    printf("Phi3(%d, %d) = %g\n", 1, 1, EvalLin2D(b, 3, 1, 1));
}
