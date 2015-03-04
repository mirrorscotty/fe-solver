#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "matrix.h"
#include "basis.h"
#include "integrate.h"
#include "mesh1d.h"
#include "isoparam.h"
#include "finite-element1d.h"

#include "material-data.h"

#include "heat-transfer.h"
#include "diffusion.h"
#include "common.h"

extern choi_okos *comp_global;

/**
 * Calculate an average diffusivity for the entire drying process. The
 * temperature is taken to be the ambient temperature, and the average between
 * the initial and final moisture content is used.
 * @param Xdb Moisture content (Not used)
 * @param T Temperature (Not used)
 * @returns Diffusivity [m^2/s]
 */
double DiffAvg(double Xdb, double T)
{
    return DiffCh10((CINIT+CAMB)/2, TAMB);
}

/**
 * This takes the function that calculates the residual and integrates it
 * to form the element-level matrix used in solving the problem. This whole
 * function is given to the finite element solver and run repeatedly to generate
 * the global matrix given to the matrix solver.
 * @param p Finite element problem structure
 * @param elem Element to calculate the matrix for
 * @param guess Estimated value of the solution at the current time step
 *
 * @returns An n*v x n*v matrix containing the integrated values of the
 *      residual, where "n" is the number of basis functions used. For a
 *      linear basis, n = 2. The value for "v" is the number of PDEs being
 *      solved simultaneously.
 */
matrix* CreateElementMatrix(struct fe1d *p, Elem1D *elem, matrix *guess)
{
    basis *b;
    b = p->b;
    
    int v = p->nvars;
    
    int i, j;
    double value = 0;
    matrix *m;
    
    m = CreateMatrix(b->n*v, b->n*v);
    
    for(i=0; i<b->n*v; i+=v) {
        for(j=0; j<b->n*v; j+=v) {
            #ifdef TVAR
            value = quad1d3generic(p, guess, elem, &ResHeat, i/v, j/v);
            setval(m, value, i+TVAR, j+TVAR);
            #endif
            #ifdef CVAR
            value = quad1d3generic(p, guess, elem, &ResMass, i/v, j/v);
            setval(m, value, i+CVAR, j+CVAR);
            #endif
            #ifdef VAP_MODEL
            value = quad1d3generic(p, guess, elem, &ResVap, i/v, j/v);
            setval(m, value, i+CVAR, j+CVAR);
            #endif
        }
    }

    return m;
}

/**
 * This takes the function that calculates the element-level coefficient matrix
 * multiplying \f$\frac{\partial T_i}{\partial t}\f$. As with the
 * CreateElementMatrix function, this is passed to the FEM solver and run
 * repeatedly to generate the global matrix.
 * @param p Finite element problem structure
 * @param elem Element to calculate the matrix for
 * @param guess Estimated value of the solution at the current time step
 *
 * @returns An n*v x n*v matrix containing the integrated values of the
 *      portion of the residual multiplying
 *      \f$\frac{\partial T_i}{\partial t}\f$, where "n" is the number of basis
 *      functions used. For a linear basis, n = 2. The value for "v" is the
 *      number of PDEs being solved simultaneously.
 * @see CreateElementMatrix
 */
matrix* CreateDTimeMatrix(struct fe1d *p, Elem1D *elem, matrix *guess) {
    basis *b;
    b = p->b;
    
    int v = p->nvars;
    
    int i, j;
    double value = 0;
    matrix *m;
    
    m = CreateMatrix(b->n*v, b->n*v);
    
    for(i=0; i<b->n*v; i+=v) {
        for(j=0; j<b->n*v; j+=v) {
            #ifdef TVAR
            value = quad1d3generic(p, guess, elem, &ResDtHeat, i/v, j/v);
            setval(m, value, i+TVAR, j+TVAR);
            #endif
            #ifdef CVAR
            value = quad1d3generic(p, guess, elem, &ResDtMass, i/v, j/v);
            setval(m, value, i+CVAR, j+CVAR);
            #endif
            #ifdef VAP_MODEL
            value = quad1d3generic(p, guess, elem, &ResDtVap, i/v, j/v);
            setval(m, value, i+CVAR, j+CVAR);
            #endif
        }
    }

    return m;
}

/**
 * Creates an integrated, element-level load vector (source term + boundary
 * conditions). The boundary conditions are added separately, and the source
 * term for this equation is equal to zero, so this simply returns a column
 * matrix of zeros of the appropriate size.
 */
matrix* CreateElementLoad(struct fe1d *p, Elem1D *elem, matrix *guess) {
    basis *b;
    b = p->b;
    
    int v = p->nvars;
    
    matrix *m;
    
    m = CreateMatrix(b->n*v, 1);

    return m;
}

/**
 * Returns true if the specified row (node) is on the right-most boundary of the
 * domain.
 * @param p Finite element problem structure
 * @param row Row to check to see if it's on the boundary. Since this is a 1D
 *      problem, the row number corresponds to the node number.
 * @returns True if the node is on the right-hand boundary and false otherwise.
 */
int IsOnRightBoundary(struct fe1d *p, int row)
{
    if(row == len(p->mesh->nodes)-1)
        return 1;
    else
        return 0;
}

/**
 * Returns true if the specified row (node) is on the left-most boundary of the
 * domain.
 * @param p Finite element problem structure
 * @param row Row to check to see if it's on the boundary. Since this is a 1D
 *      problem, the row number corresponds to the node number.
 * @returns True if the node is on the left-hand boundary and false otherwise.
 */
int IsOnLeftBoundary(struct fe1d *p, int row)
{
    if(row == 0)
        return 1;
    else
        return 0;
}

/**
 * Take the tests functions that define where the boundaries are and the
 * functions for boundary conditions and apply them to the problem. This is
 * given to the FEM solver and run each time the global matrices are
 * recalculated.
 * @param p Finite element problem structure
 */
void ApplyAllBCs(struct fe1d *p)
{
#ifdef TVAR
    double Bi = BiotNumber(p->charvals);
#endif
#ifdef CVAR
    double Bim = BiotNumber(p->chardiff);
#endif
    
    /* BC at x=L:
     * This approximates any Biot number larger than 100 as Bi->infty. This is a
     * good approximation for this problem since it results in the outside of
     * the can reaching the external temperature incredibly quickly (<1sec). */
#ifdef TVAR
    if(Bi<100.00)
        ApplyNaturalBC1D(p, TVAR, &IsOnRightBoundary, &ConvBCHeat);
    else
        ApplyEssentialBC1D(p, TVAR, &IsOnRightBoundary, &ExternalTemp);
#endif

    /* Do the same for the mass transfer boundary condition. */
#ifdef CVAR
    if(Bim<100.00)
        ApplyNaturalBC1D(p, CVAR, &IsOnRightBoundary, &ConvBCMass);
    else
        ApplyEssentialBC1D(p, CVAR, &IsOnRightBoundary, &ExternalConc);
#endif

#ifdef VAP_MODEL
    if(Bim<100.00)
        ApplyNaturalBC1D(p, PVAR, &IsOnRightBoundary, &ConvBCVap);
    else
        ApplyEssentialBC1D(p, PVAR, &IsOnRightBoundary, &ExternalConc);
#endif
    return;
}

/* Calculate the value of the deformation gradient at the specified time and
 * spatial coordinates.
 * \f[ \underline{\underline{F}}(\underline{X}, T)\f]
 * This value is calculated from the density Choi-Okos equations, and supplied
 * to the moving mesh function to recalculate the nodal values at each time
 * step.
 * @param p Finite element problem structure
 * @param X Spatial coordinate (global)
 * @param t Time step number
 * @returns Calculated value for the deformation gradient
 */
double DeformationGrad(struct fe1d *p, double X, double t)
{
    solution *s0, *sn;
    double rho0, rhon;
    double T0 = TINIT, Tn = TINIT;
    
#ifdef CVAR
    double C0, Cn;
    choi_okos *cowet0, *cowetn;
#endif
    
    s0 = FetchSolution(p, 0);
    sn = FetchSolution(p, t);

#ifdef TVAR
    Tn = uscaleTemp(p->charvals, EvalSoln1DG(p, TVAR, sn, X, 0));
    T0 = uscaleTemp(p->charvals, EvalSoln1DG(p, TVAR, s0, X, 0));
#endif
#ifdef CVAR
    Cn = uscaleTemp(p->chardiff, EvalSoln1DG(p, CVAR, sn, X, 0));
    C0 = uscaleTemp(p->chardiff, EvalSoln1DG(p, CVAR, s0, X, 0));

    cowet0 = AddDryBasis(comp_global, C0);
    cowetn = AddDryBasis(comp_global, Cn);

    rho0 = rho(cowet0, T0);
    rhon = rho(cowetn, Tn);

    DestroyChoiOkos(cowet0);
    DestroyChoiOkos(cowetn);
#else
    rhon = rho(comp_global, Tn);
    rho0 = rho(comp_global, T0);
#endif

    return rho0/rhon;
}

