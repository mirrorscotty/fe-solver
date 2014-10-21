#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "basis.h"
#include "integrate.h"
#include "mesh2d.h"
#include "isoparam.h"
#include "finite-element.h"

/* Chroniker delta function */
int ChronDelta(int a, int b)
{
    return (a==b)?1:0;
}

/* Deformation gradient tensor 
 * I'm not sure this function makes sense */
double FTensor(struct fe *p, Elem2D *elem, matrix *guess, int a, int b, double x, double y)
{
    double Fab;
    int i;
    
    Fab = ChronDelta(a, b);
    
    for(i=b; i<nRows(guess); i+=2) {
        if(a) {
            Fab += IEvalLin2Dy(p, elem, i/2, x, y) * val(guess, i, 0);
        } else {
            Fab += IEvalLin2Dx(p, elem, i/2, x, y) * val(guess, i, 0);
        }
    }
    return Fab;
}

/* Calculate the desired element of the Lagrangian strain tensor. */
double ETensor(struct fe *p, Elem2D *elem, matrix *guess, int a, int b, double x, double y)
{
    return (FTensor(p, elem, guess, b, a, x, y)*FTensor(p, elem, guess, a, b, x, y) - ChronDelta(a, b))/2;
}

/* Return the a,b component of the second Piola-Kirchhoff stress tensor */
double STensor(struct fe *p, Elem2D *elem, matrix *guess, int a, int b, double x, double y)
{
    double C = 2e3;
    double v = 0.3;
    double Sab = 0;
    Sab = C*v/((1+v)*(1-2*v)) * (ETensor(p, elem, guess, 0, 0, x, y)+ETensor(p, elem, guess, 1, 1, x, y)) * ChronDelta(a,b);
    Sab += C/(1+v)*ETensor(p, elem, guess, a, b, x, y);
    
    return Sab;
}

double ElemJdRxdx(struct fe *p, matrix *guess, Elem2D *elem, double x, double y, int f1, int f2)
{
    double value;
    double C = 2e3;
    double v = 0.3;
    double a = C*v/((1+v)*(1-2*v));
    double c = C/(1+v);

    double Fxx = FTensor(p, elem, guess, 0, 0, x, y);
    double Fxy = FTensor(p, elem, guess, 0, 1, x, y);
   // double Fyx = FTensor(p, elem, guess, 1, 0, x, y);
   // double Fyy = FTensor(p, elem, guess, 1, 1, x, y);
    
    //printf("F:\n[%g, %g; %g, %g]\n", Fxx, Fxy, Fyx, Fyy);
    
    double Sxx = STensor(p, elem, guess, 0, 0, x, y);
    double Sxy = STensor(p, elem, guess, 0, 1, x, y);
    double Syx = STensor(p, elem, guess, 1, 0, x, y);
    double Syy = STensor(p, elem, guess, 1, 1, x, y);

    //printf("S:\n[%g, %g; %g, %g]\n", Sxx, Sxy, Syx, Syy);
    
    double A = IEvalLin2Dx(p, elem, f1, x, y) * Sxx * IEvalLin2Dx(p, elem, f2, x, y);
    A += IEvalLin2Dx(p, elem, f1, x, y) * Sxy * IEvalLin2Dy(p, elem, f2, x, y);
    A += IEvalLin2Dy(p, elem, f1, x, y) * Syx * IEvalLin2Dx(p, elem, f2, x, y);
    A += IEvalLin2Dy(p, elem, f1, x, y) * Syy * IEvalLin2Dy(p, elem, f2, x, y);

    value = 0;
    value += A;
    value += a*(Fxx*IEvalLin2Dx(p, elem, f1, x, y) + Fxy*IEvalLin2Dy(p, elem, f1, x, y)) * 
               (Fxx*IEvalLin2Dx(p, elem, f2, x, y) + Fxy*IEvalLin2Dy(p, elem, f2, x, y));
    value += c/2 * (Fxx*Fxx + Fxy*Fxy)*(IEvalLin2Dx(p, elem, f1, x, y) * IEvalLin2Dx(p, elem, f2, x, y) +
                IEvalLin2Dy(p, elem, f1, x, y) * IEvalLin2Dy(p, elem, f2, x, y));
    value += c/2 *(Fxx*IEvalLin2Dx(p, elem, f2, x, y) + Fxy*IEvalLin2Dy(p, elem, f2, x, y)) * 
                  (Fxx*IEvalLin2Dx(p, elem, f1, x, y) + Fxy*IEvalLin2Dy(p, elem, f1, x, y));
    
    value *= IMapJ(p, elem, x, y);
    
    return value;
}

double ElemJdRxdy(struct fe *p, matrix *guess, Elem2D *elem, double x, double y, int f1, int f2)
{
    double value;
    double C = 2e3;
    double v = 0.3;
    double a = C*v/((1+v)*(1-2*v));
    double c = C/(1+v);

    double Fxx = FTensor(p, elem, guess, 0, 0, x, y);
    double Fxy = FTensor(p, elem, guess, 0, 1, x, y);
    double Fyx = FTensor(p, elem, guess, 1, 0, x, y);
    double Fyy = FTensor(p, elem, guess, 1, 1, x, y);
    
    value = 0;
    value += a*(Fxx*IEvalLin2Dx(p, elem, f1, x, y) + Fxy*IEvalLin2Dy(p, elem, f1, x, y)) * 
               (Fyx*IEvalLin2Dx(p, elem, f2, x, y) + Fyy*IEvalLin2Dy(p, elem, f2, x, y));
    value += c/2 * (Fyx*Fxx+Fyx*Fxy)*(IEvalLin2Dx(p, elem, f1, x, y) * IEvalLin2Dx(p, elem, f2, x, y) +
                            IEvalLin2Dy(p, elem, f1, x, y) * IEvalLin2Dy(p, elem, f2, x, y));
    value += c/2 *(Fxx*IEvalLin2Dx(p, elem, f2, x, y) + Fxy*IEvalLin2Dy(p, elem, f2, x, y)) * 
                  (Fyx*IEvalLin2Dx(p, elem, f1, x, y) + Fyy*IEvalLin2Dy(p, elem, f1, x, y));
                  
    value *= IMapJ(p, elem, x, y);

    return value;
}

double ElemJdRydx(struct fe *p, matrix *guess, Elem2D *elem, double x, double y, int f1, int f2)
{
    double value;
    double C = 2e3;
    double v = 0.3;
    double a = C*v/((1+v)*(1-2*v));
    double c = C/(1+v);

    double Fxx = FTensor(p, elem, guess, 0, 0, x, y);
    double Fxy = FTensor(p, elem, guess, 0, 1, x, y);
    double Fyx = FTensor(p, elem, guess, 1, 0, x, y);
    double Fyy = FTensor(p, elem, guess, 1, 1, x, y);
    
    value = 0;
    value += a*(Fyx*IEvalLin2Dx(p, elem, f1, x, y) + Fyy*IEvalLin2Dy(p, elem, f1, x, y)) * 
               (Fxx*IEvalLin2Dx(p, elem, f2, x, y) + Fxy*IEvalLin2Dy(p, elem, f2, x, y));
    value += c/2 * (Fxx*Fyx+Fxy*Fyy)*(IEvalLin2Dx(p, elem, f1, x, y) * IEvalLin2Dx(p, elem, f2, x, y) +
                            IEvalLin2Dy(p, elem, f1, x, y) * IEvalLin2Dy(p, elem, f2, x, y));
    value += c/2 *(Fyx*IEvalLin2Dx(p, elem, f2, x, y) + Fyy*IEvalLin2Dy(p, elem, f2, x, y)) * 
                  (Fxx*IEvalLin2Dx(p, elem, f1, x, y) + Fxy*IEvalLin2Dy(p, elem, f1, x, y));
 
    value *= IMapJ(p, elem, x, y);
    
    return value;
}

double ElemJdRydy(struct fe *p, matrix *guess, Elem2D *elem, double x, double y, int f1, int f2)
{
    double value;
    double C = 2e3;
    double v = 0.3;
    double a = C*v/((1+v)*(1-2*v));
    double c = C/(1+v);

    //double Fxx = FTensor(p, elem, guess, 0, 0, x, y);
    //double Fxy = FTensor(p, elem, guess, 0, 1, x, y);
    double Fyx = FTensor(p, elem, guess, 1, 0, x, y);
    double Fyy = FTensor(p, elem, guess, 1, 1, x, y);
    
    double Sxx = STensor(p, elem, guess, 0, 0, x, y);
    double Sxy = STensor(p, elem, guess, 0, 1, x, y);
    double Syx = STensor(p, elem, guess, 1, 0, x, y);
    double Syy = STensor(p, elem, guess, 1, 1, x, y);

    double A = IEvalLin2Dx(p, elem, f1, x, y) * Sxx * IEvalLin2Dx(p, elem, f2, x, y);
    A += IEvalLin2Dx(p, elem, f1, x, y) * Sxy * IEvalLin2Dy(p, elem, f2, x, y);
    A += IEvalLin2Dy(p, elem, f1, x, y) * Syx * IEvalLin2Dx(p, elem, f2, x, y);
    A += IEvalLin2Dy(p, elem, f1, x, y) * Syy * IEvalLin2Dy(p, elem, f2, x, y);

    value = 0;
    value += A;
    value += a*(Fyx*IEvalLin2Dx(p, elem, f1, x, y) + Fyy*IEvalLin2Dy(p, elem, f1, x, y)) * 
               (Fyx*IEvalLin2Dx(p, elem, f2, x, y) + Fyy*IEvalLin2Dy(p, elem, f1, x, y));
    value += c/2 * (Fyx*Fyx+Fyy*Fyy)*(IEvalLin2Dx(p, elem, f1, x, y) * IEvalLin2Dx(p, elem, f2, x, y) +
                            IEvalLin2Dy(p, elem, f1, x, y) * IEvalLin2Dy(p, elem, f2, x, y));
    value += c/2 *(Fyx*IEvalLin2Dx(p, elem, f2, x, y) + Fyy*IEvalLin2Dy(p, elem, f2, x, y)) * 
                  (Fyx*IEvalLin2Dx(p, elem, f1, x, y) + Fyy*IEvalLin2Dy(p, elem, f1, x, y));

    value *= IMapJ(p, elem, x, y);
    
    return value;
}

matrix* CreateElementMatrix(struct fe *p, Elem2D *elem, matrix *guess)
{
    basis *b;
    b = p->b;
    
    int nvars = p->nvars;
    
    int i, j;
    double value = 0;
    matrix *m;
    
    m = CreateMatrix(b->n*nvars, b->n*nvars);
    
    for(i=0; i<b->n*nvars; i+=nvars) {
        for(j=0; j<b->n*nvars; j+=nvars) {
            /* dRx/dx */
            value = quad2d3generic(p, guess, elem, &ElemJdRxdx, i/2, j/2);
            setval(m, value, i, j);
            
            /* dRy/dx */
            value = quad2d3generic(p, guess, elem, &ElemJdRydx, i/2, j/2);
            setval(m, value, i+1, j);
            
            /* dRx/dy */
            value = quad2d3generic(p, guess, elem, &ElemJdRxdy, i/2, j/2);
            setval(m, value, i, j+1);
            
            /* dRy/dy */
            value = quad2d3generic(p, guess, elem, &ElemJdRydy, i/2, j/2);
            setval(m, value, i+1, j+1);
        }
    }

    return m;
}

/* Create the load vector */
matrix* CreateElementLoad(struct fe *p, Elem2D *elem, matrix *guess) {
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

int IsOnYAxis(struct fe *p, int row)
{
    double x = valV(p->mesh->nodes[row/p->nvars], 0);

    if(fabs(x) < 1e-5)
        return 1;
    else
        return 0;
}

int IsOnXAxis(struct fe *p, int row)
{
    double y = valV(p->mesh->nodes[row/p->nvars], 1);

    if(fabs(y) < 1e-5)
        return 1;
    else
        return 0;
}

int IsOnRightBoundary(struct fe *p, int row)
{
    double width = 1;
    double x = valV(p->mesh->nodes[row/p->nvars], 0);
    
    if(fabs(x - width) < 1e-5)
        return 1;
    else
        return 0;
}

int IsOnTopBoundary(struct fe *p, int row)
{
    double height = p->mesh->y2 - p->mesh->y1;
    double y = valV(p->mesh->nodes[row/p->nvars], 1);
    
    if(fabs(y - height) < 1e-5)
        return 1;
    else
        return 0;
}

double Zero(struct fe *p, int row)
{
    return 0.0;
}

double TopPressure(struct fe* p, int row)
{
    double dx = p->mesh->x2 - p->mesh->x1;
    dx = dx/p->mesh->nelemx;
    double Pressure = p->P;
    return Pressure * dx/2;
}

double LeftPressure(struct fe* p, int row)
{
    double dy = p->mesh->y2 - p->mesh->y1;
    double Pressure = p->P;
    double alpha = p->a;
    return Pressure * alpha * dy/2;
}

void ApplyAllBCs(struct fe *p)
{
    
    // BC at x=0
    ApplyNaturalBC(p, 0, &IsOnYAxis, &LeftPressure);
    
    // BC at y=H
    ApplyNaturalBC(p, 1, &IsOnTopBoundary, &TopPressure);
    
    // BC at y=0
    ApplyEssentialBC(p, 0, &IsOnXAxis, &Zero);
    ApplyEssentialBC(p, 1, &IsOnXAxis, &Zero);
    return;
}

/* Functions to format the results */
matrix* GetXDisplacements(matrix *Soln) {
    int i;
    matrix* xd;
    xd = CreateMatrix(nRows(Soln), 1);
    
    for(i=0; i<nRows(xd); i++) {
        setval(xd, val(Soln, i*2, 0), i, 0);
    }
    
    return xd;
}

matrix* GetYDisplacements(matrix *Soln)
{
    int i;
    matrix* yd;
    yd = CreateMatrix(nRows(Soln), 1);
    
    for(i=0; i<nRows(yd); i++) {
        setval(yd, val(Soln, i*2+1, 0), i, 0);
    }
    
    return yd;
}

matrix* GetDeformedCoords(struct fe *p, matrix *Soln)
{
    matrix *Def;
    int i;
    Def = CreateMatrix(nRows(Soln)/2, 2);
    
    for(i=0; i<nRows(Soln)/2; i++) {
        setval(Def, val(Soln, 2*i, 0) + valV(GetNodeCoordinates(p->mesh, i), 0),  i, 0);
        setval(Def, val(Soln, 2*i+1, 0) + valV(GetNodeCoordinates(p->mesh, i), 1), i, 1);
    }
    
    return Def;
}

matrix* FormatDisplacements(struct fe *p, matrix *Soln)
{
    matrix *Def;
    int i;
    Def = CreateMatrix(nRows(Soln)/2, 2);
    
    for(i=0; i<nRows(Soln)/2; i++) {
        setval(Def, val(Soln, 2*i, 0),  i, 0);
        setval(Def, val(Soln, 2*i+1, 0), i, 1);
    }
    
    return Def;
}
    
int main(int argc, char *argv[])
{
    Mesh2D *mesh;
    basis *b;
    matrix *E;
    struct fe* problem;

    /* Make a linear 2D basis */
    b = MakeLinBasis(2);

    /* Create a uniform mesh */
    mesh = GenerateUniformMesh2D(0.0, 1.0,
                                 0.0, 2.0,
                                 8, 8);
    
    problem = CreateFE(b, mesh, &CreateElementMatrix, &CreateElementLoad, &ApplyAllBCs);
    problem->nvars = 2;
    
    problem->P = .1;
    problem->a = .01;
    
    E = NLinSolve(problem, NULL);

    //mtxprnt(E);
    //puts("");
    mtxprnt(GetDeformedCoords(problem, E));
    
    return 0;
}
