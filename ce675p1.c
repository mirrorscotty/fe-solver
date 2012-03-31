#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "basis.h"
#include "mtxsolver.h"
#include "integrate.h"
#include "mesh.h"
#include "finite-element.h"

int ChronDelta(int a, int b)
{
    return (a==b)?1:0;
}
double FTensor(struct fe *p, matrix *guess, int a, int b)
{
    double Fab = 0;
    int i;
    
    Fab = ChronDelta(a, b);
    
    for(i=b; i<mtxlen2(guess); i+=2) {
        if(a) {
            Fab += quad2d3(p->b, i, 0, 1) * val(guess, i, 0);
        } else {
            Fab += quad2d3(p->b, i, 1, 0) * val(guess, i, 0);
        }
    }
    return Fab;
}

double ETensor(struct fe *p, matrix *guess, int a, int b)
{
    return (FTensor(p, guess, b, a)*FTensor(p, guess, a, b) - ChronDelta(a, b))/2;
}

double STensor(struct fe *p, matrix *guess, int a, int b)
{
    double C = 2e3;
    double v = 0.3;
    double Sab = 0;
    Sab = C*v/((1+v)*(1-2*v)) * (ETensor(p, guess, 0, 0)+ETensor(p, guess, 1, 1)) * ChronDelta(a,b);
    Sab += C/(1+v)*ETensor(p, guess, a, b);
    
    return Sab;
}

matrix* CreateElementMatrix(struct fe *p, matrix *guess)
{
    basis *b;
    Elem2D *elem;
    
    b = p->b;
    int nvars = p->nvars;
    int var;
    
    /* Todo: Fix this. */
    elem = &(p->mesh->elem[0]);
    
    int i, j;
    double value = 0;
    double dx = elem->dx;
    double dy = elem->dy;
    matrix *m;
    
    
    double Fxx = FTensor(p, guess, 0, 0);
    double Fxy = FTensor(p, guess, 0, 1);
    double Fyx = FTensor(p, guess, 1, 0);
    double Fyy = FTensor(p, guess, 1, 1);
    
    //printf("F:\n%g %g\n%g %g\n\n", Fxx, Fxy, Fyx, Fyy);
    
    double Sxx = STensor(p, guess, 0, 0);
    double Sxy = STensor(p, guess, 0, 1);
    double Syx = STensor(p, guess, 1, 0);
    double Syy = STensor(p, guess, 1, 1);

    //printf("S:\n%g %g\n%g %g\n\n", Sxx, Sxy, Syx, Syy);
    
    m = CreateMatrix(b->n*nvars, b->n*nvars);
    
    double C = 2e3;
    double v = 0.3;
    double a = C*v/((1+v)*(1-2*v));
    double c = C/(1+v);
    
    for(i=0; i<b->n*nvars; i+=nvars) {
        for(j=0; j<b->n*nvars; j+=nvars) {
            /* Scalars */
            double A = quad2d3(b, i/2, 1, 0) * Sxx * quad2d3(b, j/2, 1, 0);
            A += quad2d3(b, i/2, 1, 0) * Sxy * quad2d3(b, j/2, 0, 1);
            A += quad2d3(b, i/2, 0, 1) * Syx * quad2d3(b, j/2, 1, 0);
            A += quad2d3(b, i/2, 0, 1) * Syy * quad2d3(b, j/2, 0, 1);
            
            /* dRx/dx */
            value = 0;
            value += A;
            value += a*(Fxx*quad2d3(b, i/2, 1, 0) + Fxy*quad2d3(b, i/2, 0, 1)) * 
                       (Fxx*quad2d3(b, j/2, 1, 0) + Fxy*quad2d3(b, i/2, 0, 1));
            value += c/2 * Fxx*Fxx*(quad2d3(b, i/2, 1, 0) * quad2d3(b, j/2, 1, 0) +
                                    quad2d3(b, i/2, 0, 1) * quad2d3(b, j/2, 0, 1));
            value += c/2 *(Fxx*quad2d3(b, i/2, 1, 0) + Fxy*quad2d3(b, i/2, 0, 1)) * 
                          (Fxx*quad2d3(b, j/2, 1, 0) + Fxy*quad2d3(b, i/2, 0, 1));
            setval(m, value, i, j);
            
            /* dRy/dx */
            value = 0;
            value += a*(Fyx*quad2d3(b, i/2, 1, 0) + Fyy*quad2d3(b, i/2, 0, 1)) * 
                       (Fxx*quad2d3(b, j/2, 1, 0) + Fxy*quad2d3(b, i/2, 0, 1));
            value += c/2 * Fyx*Fxy*(quad2d3(b, i/2, 1, 0) * quad2d3(b, j/2, 1, 0) +
                                    quad2d3(b, i/2, 0, 1) * quad2d3(b, j/2, 0, 1));
            value += c/2 *(Fxx*quad2d3(b, i/2, 1, 0) + Fxy*quad2d3(b, i/2, 0, 1)) * 
                          (Fxx*quad2d3(b, j/2, 1, 0) + Fxy*quad2d3(b, i/2, 0, 1));
            setval(m, value, i+1, j);
            
            /* dRx/dy */
            value = 0;
            value += a*(Fxx*quad2d3(b, i/2, 1, 0) + Fxy*quad2d3(b, i/2, 0, 1)) * 
                       (Fyx*quad2d3(b, j/2, 1, 0) + Fyy*quad2d3(b, i/2, 0, 1));
            value += c/2 * Fxy*Fyx*(quad2d3(b, i/2, 1, 0) * quad2d3(b, j/2, 1, 0) +
                                    quad2d3(b, i/2, 0, 1) * quad2d3(b, j/2, 0, 1));
            value += c/2 *(Fxx*quad2d3(b, i/2, 1, 0) + Fxy*quad2d3(b, i/2, 0, 1)) * 
                          (Fxx*quad2d3(b, j/2, 1, 0) + Fxy*quad2d3(b, i/2, 0, 1));
            setval(m, value, i, j+1);
            
            /* dRy/dy */
            value = 0;
            value += A;
            value += a*(Fyx*quad2d3(b, i/2, 1, 0) + Fyy*quad2d3(b, i/2, 0, 1)) * 
                       (Fyx*quad2d3(b, j/2, 1, 0) + Fyy*quad2d3(b, i/2, 0, 1));
            value += c/2 * Fyy*Fyy*(quad2d3(b, i/2, 1, 0) * quad2d3(b, j/2, 1, 0) +
                                    quad2d3(b, i/2, 0, 1) * quad2d3(b, j/2, 0, 1));
            value += c/2 *(Fxx*quad2d3(b, i/2, 1, 0) + Fxy*quad2d3(b, i/2, 0, 1)) * 
                          (Fxx*quad2d3(b, j/2, 1, 0) + Fxy*quad2d3(b, i/2, 0, 1));
            setval(m, value, i+1, j+1);
        }
    }

    return m;
}

/* Create the load vector */
matrix* CreateElementLoad(struct fe *p, matrix *guess) {
    int i;
    matrix *f;
    
    basis *b;
    b = p->b;

    f = CreateMatrix(b->n, 1);

    for(i=0; i<b->n; i++) {
        setval(f, 0, i, 0);
    }

    return f;
}

/* Simple function to alter J and F so that Dirichlet boundary conditions are
 * imposed on both ends of the domain. In this case, "leftbc" is imposed at
 * c = 0, and rightbc is imposed at c = 1.
 */
void ApplyDBC(matrix* J, matrix* F, basis *b, int node, double value)
{
    int rows = mtxlen2(J);
    int i;

    for(i=0;i<rows; i++) {
        setval(J, 0, node, i);
        //setval(J, 0, rows-o, i);
    }

    setval(J, 1, node, node);

    setval(F, value, node, 0);
}

void ApplyNBC(struct fe *problem, int node, double Pressure)
{
    double P;
    if(node%2) {
        P = Pressure/problem->mesh->elem[0].dy/2;
    } else {
        P = Pressure/problem->mesh->elem[0].dx/2;
    }
    
    addval(problem->F, P, node, 0);
}

void ApplyAllBCs(struct fe *p)
{
    double P = p->P;
    double a = p->a;
    matrix *J, *F;
    basis *b;
    Mesh2D *mesh;
    
    b = p->b;
    mesh = p->mesh;
    F = p->F;
    J = p->J;
    
    int i;
    
    // BC at y=0
    for(i=0; i< (mesh->nelemx+1)*(mesh->nelemy+1)*p->nvars; i+= p->nvars*(mesh->nelemy+1)) {
        ApplyDBC(J, F, b, i, 0 );
        ApplyDBC(J, F, b, i+1, 0);
        //ApplyDBC(J, F, b, 4, 0);
        //ApplyDBC(J, F, b, 5, 0);
    }
    
    // BC at x=0
    for(i=p->nvars; i< (mesh->nelemy+1)*p->nvars; i+=p->nvars) {
        //ApplyNBC(p, 2, P*a);
        ApplyNBC(p, i, P*a);
    }
    
    // BC at y=H
    for(i=(mesh->nelemy+1)*p->nvars; i< (mesh->nelemx+1)*(mesh->nelemy+1)*p->nvars; i+=p->nvars*(mesh->nelemy+1)) {
//        ApplyNBC(p, 3, P);
//        ApplyNBC(p, 7, P);
        ApplyNBC(p, i+1, P);
    }
        
        
        
        
 //   }
    return;
}

/* Functions to format the results */
matrix* GetXDisplacements(matrix *Soln) {
    int i;
    matrix* xd;
    xd = CreateMatrix(mtxlen2(Soln), 1);
    
    for(i=0; i<mtxlen2(xd); i++) {
        setval(xd, val(Soln, i*2, 0), i, 0);
    }
    
    return xd;
}

matrix* GetYDisplacements(matrix *Soln)
{
    int i;
    matrix* yd;
    yd = CreateMatrix(mtxlen2(Soln), 1);
    
    for(i=0; i<mtxlen2(yd); i++) {
        setval(yd, val(Soln, i*2+1, 0), i, 0);
    }
    
    return yd;
}

matrix* GetDeformedCoords(struct fe *p, matrix *Soln)
{
    matrix *Def;
    int i;
    Def = CreateMatrix(mtxlen2(Soln)/2, 2);
    
    for(i=0; i<mtxlen2(Soln)/2; i++) {
        setval(Def, val(Soln, 2*i, 0) + valV(GetNodeCoordinates(p->mesh, i), 0), i, 0);
        setval(Def, val(Soln, 2*i+1, 0) + valV(GetNodeCoordinates(p->mesh, i), 1), i, 1);
    }
    
    return Def;
}
    
int main(int argc, char *argv[])
{
    Mesh2D *mesh;
    basis *b;

    matrix *J, *F, *E;
    
    struct fe* problem;

    b = MakeLinBasis(2);

    /* Create a uniform mesh */
    mesh = GenerateUniformMesh2D(0.0, 1.0,
                                 0.0, 1.0,
                                 2, 2);
    
    problem = CreateFE(b, mesh, &CreateElementMatrix, &CreateElementLoad, &ApplyAllBCs);
    problem->nvars = 2;
    
    problem->P = 1;
    problem->a = 0;
    
    //AssembleJ(problem, NULL);
    
    //mtxprnt(problem->J);
    
    //E = SolveMatrixEquation(problem->J, problem->F);
    //E = NLinSolve(problem, NULL);
    //AssembleJ(problem, NULL);
    //problem->F = CreateMatrix(mtxlen2(problem->J), 1);
    
    //mtxprnt(problem->J);
    
    E = NLinSolve(problem, NULL);

    mtxprnt(E);
    puts("");
    mtxprnt(GetDeformedCoords(problem, E));
    
    

    return 0;
}
