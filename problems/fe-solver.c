#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "basis.h"
#include "integrate.h"

/* Equation: u'' = Pe*u' */
/* Create the element matrix for the specified values of Pe and h */
matrix* CreateElementMatrixOld(double Pe, double h)
{
    matrix *elem;
 
    elem = CreateMatrix(2, 2);
    setval(elem, 1/h-Pe/2, 0, 0);
    setval(elem, -1/h+Pe/2, 0, 1);
    setval(elem, -1/h-Pe/2, 1, 0);
    setval(elem, 1/h+Pe/2, 1, 1);

    return elem;
}

matrix* CreateElementMatrix(basis *b, double Pe, double h)
{
    int i, j;
    double value;
    matrix *elem;

    elem = CreateMatrix(b->n, b->n);

    for(i=0; i<b->n; i++) {
        for(j=0; j<b->n; j++) {
            value = -1/h*quad300(b->dphi[i], b->dphi[j]);
            value += Pe*quad300(b->phi[j], b->dphi[i]);
            setval(elem, value, i, j);
            if((b->n == 4) && (i%2 == 1)) {
                setval(elem, val(elem, i, j)*h, i, j);
            }
        }
    }

    return elem;
}

/* Create the load vector */
matrix* CreateElementLoad(basis* b, double Pe, double h) {
    int i;
    matrix *f;

    f = CreateMatrix(b->n, 1);

    for(i=0; i<b->n; i++) {
        setval(f, 0, i, 0);
    }

    return f;
}

/* Equation: u'' + u = 1 */
/* Test functions to create the element and load matricies based on the equation
 * in Problem 1
 */
matrix* testelem(basis* b, double Pe, double h)
{
    matrix *elem;
    int i, j;
 
    elem = CreateMatrix(b->n, b->n);
    for(i=0; i<b->n; i++) {
        for(j=0; j<b->n; j++) {
            setval(elem, -1/h*quad300(b->dphi[i], b->dphi[j]) + h*quad300(b->phi[i], b->phi[j]), i, j);
            if((b->n == 4) && (i%2 == 1)) {
                setval(elem, val(elem, i, j)*h, i, j);
            }
        }
    }
    return elem;
}

matrix* testload(basis* b, double Pe, double h) {
    matrix *f;
    int i;

    f = CreateMatrix(b->n, 1);

    for(i=0; i<b->n; i++) {
        setval(f, h*quad3(b->phi[i]), i, 0);
        if((b->n == 4) && (i%2 == 1)) {
            setval(f, val(f, i, 0)*h, i, 0);
        }
    }

    return f;
}

/* Assemble the global coefficient matrix for a non-uniform mesh. The first
 * argument is a function pointer to the function which creates the element
 * matrix, and the second is the mesh to be used.
 */
matrix* AssembleJ(matrix* (*makej)(basis*, double, double), basis *b, matrix* mesh, double Pe)
{
    matrix *J, *j;
    int n, i, x, y, rows;
    int r = b->n-b->overlap; /* Number of rows to add */

    /* Determine the number of elements from the mesh */
    n = nRows(mesh);
    /* Determine the number of rows required for the global matrix */
    rows = n*r+b->overlap;

    /* Create a blank global matrix */
    J = CreateMatrix(rows, rows);

    for(i=0; i<rows-r; i=i+r) {
	/* Generate the element matrix for the specified element width */
	j = makej(b, Pe, val(mesh, i/r, 0));

	/* Add the values of the element matrix to the global matrix */
        for(x=0; x<b->n; x++) {
            for(y=0; y<b->n; y++) {
                addval(J, val(j, x, y), i+x, i+y);
            }
        }

        /* Clean up */
        DestroyMatrix(j);
    }

    return J;
}

matrix* AssembleF(matrix* (*makef)(basis*, double, double), basis *b, matrix* mesh, double Pe)
{
    matrix *F, *f;
    int n, i, j, rows;
    int r = b->n - b->overlap;

    n = nRows(mesh);
    /* Determine the number of rows required for the global matrix */
    rows = n*r+b->overlap;

    F = CreateMatrix(rows, 1);

    for(i=0; i<rows-r; i=i+r) {
        f = makef(b, Pe, val(mesh, i/r, 0));

        for(j=0; j<b->n; j++) {
            addval(F, val(f, j, 0), i+j, 0);
        }

        DestroyMatrix(f);
    }

    return F;
}


/* Simple function to alter J and F so that Dirichlet boundary conditions are
 * imposed on both ends of the domain. In this case, "leftbc" is imposed at
 * c = 0, and rightbc is imposed at c = 1.
 */
void ApplyBoundaryConditions(matrix* J, matrix* F, basis *b, double leftbc, double rightbc)
{
    int rows = nRows(J);
    int i;
    int o = b->overlap;

    for(i=0;i<rows; i++) {
        setval(J, 0, 0, i);
        setval(J, 0, rows-o, i);
    }

    setval(J, 1, 0, 0);
    setval(J, 1, rows-o, rows-o);

    setval(F, leftbc, 0, 0);
    setval(F, rightbc, rows-o, 0);
}

/* Determine the element width based on the size of the domain and the number
 * of elements.
 */
matrix* GenerateUniformMesh(double left, double right, int n)
{
    double h;
    int i;
    matrix *mesh;
    mesh = CreateMatrix(n, 1);

    h = (right-left)/n;

    for(i=0; i<n; i++)
	setval(mesh, h, i, 0);

    return mesh;
}

matrix* MeshXCoords(matrix *mesh, basis* b, double left, double right)
{
    matrix *x;
    int i;
    /* Number of elements */
    int n = nRows(mesh);
    /* Nodes per element */
    int nodes = b->n-b->overlap;

    switch(b->n) {
        case 2:
        case 4:
            x = CreateMatrix(n+1, 1);
            setval(x, left, 0, 0);
            for(i=1; i<=n; i++) {
                setval(x, val(x, i-1, 0) + val(mesh, i-1, 0), i, 0);
            }
        break;
        case 3:
            x = CreateMatrix(2*n+1, 1);
            setval(x, left, 0, 0);
            for(i=1; i<=2*n; i++) {
                setval(x, val(x, i-1, 0) + val(mesh, (i-1)/2, 0)/2, i, 0);
            }
            break;
    }


    if((val(x, n, 0) - right) > 1e-5) {
        fprintf(stderr, "Warning: mesh not aligned.\n");
    }

    return x;
}

double EvalSolution(matrix *m, basis *b, double x, double left, double right, matrix *soln)
{
    int i;
    double c;
    matrix *xcoord;
    xcoord = MeshXCoords(m, b, left, right);
    switch(b->n) {
	case 2:
	    for(i=0; x>val(xcoord, i, 0)+val(m,i,0)&&i+2<nRows(xcoord); i++);
	    c = (x-val(xcoord, i, 0))/val(m, i, 0);
	    return val(soln, i, 0) * b->phi[0](c)
	            + val(soln, i+1, 0) * b->phi[1](c);
	    break;
	case 3:
	    for(i=0; x>val(xcoord, i, 0)+val(m, i/2, 0) && i+3<nRows(xcoord); i+=2);
	    c = (x-val(xcoord, i, 0))/val(m, i/2, 0);
	    return val(soln, i, 0) * b->phi[0](c) 
                + val(soln, i+1, 0) * b->phi[1](c)
                + val(soln, i+2, 0) * b->phi[2](c);
	    break;
	case 4:
	    for(i=0; x>val(xcoord, i, 0)+val(m,i,0)&&i+2<nRows(xcoord); i++);
	    c = (x-val(xcoord, i, 0))/val(m, i, 0);
	    return val(soln, 2*i, 0) * b->phi[0](c)
	            + val(soln, 2*i+1, 0) * b->phi[1](c)// * val(m, i, 0)
                + val(soln, 2*i+2, 0) * b->phi[2](c)
                + val(soln, 2*i+3, 0) * b->phi[3](c);// * val(m, i, 0);
	    break;
    }
    return 0;
}

double Exact(double x, double Pe)
{
    return (exp(Pe*x)-1)/(exp(Pe)-1);
//    return (cos(x)+1);
}

double EValue(matrix *m, basis *b, double left, double right, matrix *soln, double Pe)
{
    double x[] = {-0.774596669241483,
                   0.,
                   0.774596669241483};
    double w[] = {0.555555555555555,
                  0.888888888888888,
                  0.555555555555555}; 
    double result=0;
    int i;

    for(i=0; i<3; i++) {
	result += pow(EvalSolution(m, b, 0.5*(x[i]+1), left, right, soln) - Exact(0.5*(x[i]+1), Pe), 2) * w[i]/2;
    }

    return pow(result, .5);
}

#define PI 3.141592654

/* This program takes up to three arguments. The first is the Peclet number to
 * use in the calculation. The second is the number of elements to use when
 * solving the problem, and the third (optional) argument is the filename to
 * save the solution to. */
int main(int argc, char *argv[])
{
    double Pe, h;
    double left = 0;
    //double right = 1;
    double right = PI/2;
    int n, i;
    matrix *J, *F, *A, *mesh, *x, *E, *x1, *y, *xy;
    basis *b;

    b = MakeCubicBasis();

    /* Parse arguments */
    if(argc < 2) {
        fprintf(stderr, "Too few arguments: exiting.\n");
        return 1;
    }
    Pe = atof(argv[1]);

    n = atoi(argv[2]);

    /* Create a uniform mesh */
    mesh = GenerateUniformMesh(left, right, n);

    /* Use a mesh with these node spacings */
    //mesh = ParseMatrix("[.75;.15;.0125;.0125;.0125;.0125;.0125;.0125;.0125;.0125]");

    x = MeshXCoords(mesh, b, left, right);

    /* All the nodes are equally spaced. */
    h = val(mesh, 0, 0);

    //J = AssembleJ(&CreateElementMatrix, b, mesh, Pe);
    J = AssembleJ(&testelem, b, mesh, Pe);

    //F = AssembleF(&CreateElementLoad, b, mesh, Pe);
    F = AssembleF(&testload, b, mesh, Pe);

    //ApplyBoundaryConditions(J, F, b, 0, 1);
    ApplyBoundaryConditions(J, F, b, 2, 1);
    
    E = SolveMatrixEquation(J, F);
    printf("E = %1.10f\n", EValue(mesh, b, left, right, E, Pe));
    printf("f(pi/4) = %f\n", EvalSolution(mesh, b, PI/4, left, right, E));

    x1 = linspace(0, PI/2, 100);
    x1 = mtxtrn(x1);
    y = CreateMatrix(100, 1);
    for(i=0; i<100; i++) {
        setval(y, EvalSolution(mesh, b, val(x1, i, 0), left, right, E), i, 0);
    }
    xy = AugmentMatrix(x1, y);
    
    /* Print out the result */
    if(argc == 4) {
        mtxprntfile(xy, argv[3]);
    } else {
        mtxprnt(x);
        puts("");
        mtxprnt(E);
    }
    
    /* Clean up the allocated memory */
    DestroyMatrix(J);
    DestroyMatrix(F);
    DestroyMatrix(E);
    DestroyMatrix(mesh);

    return 0;
}
