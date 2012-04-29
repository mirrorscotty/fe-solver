/*
 * Series of test functions for figuring stuff out. Pretty much nothing here is
 * used for actual analysis.
 */

#include <stdio.h>
#include "basis.h"
#include "mesh.h"
#include "isoparam.h"
#include "integrate.h"
#include "finite-element.h"

/*
int main(int argc, char *argv[])
{
    Mesh2D *mesh;
    struct fe *p;
    int i = 5;

    p = CreateFE(MakeQuadBasis(2), mesh, NULL, NULL, NULL); 

    //mesh = GenerateUniformMesh2D(0, 10, 0, 10, 10, 10);
    mesh = MakeSpheroidMesh(p->b, 0, 5, 5, 5);
    PrintVector(mesh->elem[i]->points[0]);
    PrintVector(mesh->elem[i]->points[1]);
    PrintVector(mesh->elem[i]->points[2]);
    PrintVector(mesh->elem[i]->points[3]);

    printf("dx/dxi = %g\n", IMapXXi(p, mesh->elem[i], .5, .5));
    printf("dx/deta = %g\n", IMapXEta(p, mesh->elem[i], .5, .5));
    printf("dy/dxi = %g\n", IMapYXi(p, mesh->elem[i], .5, .5));
    printf("dy/deta = %g\n", IMapYEta(p, mesh->elem[i], .5, .5));
    MeshPrint(mesh);

    DestroyMesh2D(mesh);
    DestroyFE(p);

    return 0;
}
*/

/*
int main(int argc, char *argv[])
{
    basis *b;
    b = MakeQuadBasis(2);
    PrintQuadLocations(b);
    DestroyBasis(b);
    return 0;
}
*/

int main()
{
    Mesh2D *mesh;
    struct fe *p;
    mesh = GenerateUniformMesh2D(0, 1, 0, 1,
                                 1, 1);
    p = CreateFE(MakeLinBasis(2), mesh, NULL, NULL, NULL);
    printf("%g\n", quad2d32d3(p, mesh->elem[0], 0, 0, 0, 1));
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

/* Same thing as above for quadratic elements */
void PrintQuadLocations(basis *b)
{
    printf("Phi0(%g, %g) = %g\n", .0, .0, EvalQuad2D(b, 0, 0, 0));
    printf("Phi1(%g, %g) = %g\n", .0, .5, EvalQuad2D(b, 1, 0, .5));
    printf("Phi2(%g, %g) = %g\n", .0, 1., EvalQuad2D(b, 2, 0, 1));
    printf("Phi3(%g, %g) = %g\n", .5, 0., EvalQuad2D(b, 3, .5, 0));
    printf("Phi4(%g, %g) = %g\n", .5, .5, EvalQuad2D(b, 4, .5, .5));
    printf("Phi5(%g, %g) = %g\n", .5, 1., EvalQuad2D(b, 5, .5, 1));
    printf("Phi6(%g, %g) = %g\n", 1., 0., EvalQuad2D(b, 6, 1, 0));
    printf("Phi7(%g, %g) = %g\n", 1., .5, EvalQuad2D(b, 7, 1, .5));
    printf("Phi8(%g, %g) = %g\n", 1., 1., EvalQuad2D(b, 8, 1, 1));
}
