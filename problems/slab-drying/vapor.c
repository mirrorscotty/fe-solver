/**
 * @file diffusion.c
 * Solves the nonlinear diffusion equation.
 * \f[
 * \frac{\partial C}{\partial t} = \nabla\cdot\left(D\nabla C\right)
 * \f]
 *
 * \f[
 * \underline{\underline{F}}(\underline{X}, t) = \frac{\rho(0)}{\rho(t)}
 * \f]
 *
 * Finite Element Formulation (Global Coordinates):
 * \f[
 * R_C^i = \int_0^L \left\{
 *      \frac{\partial C_i}{\partial t}\phi_i\phi_j
 *      + C_i \frac{\partial\phi_i}{\partial t}\phi_j
 *      -C_i\frac{\partial D}{\partial x}\frac{\partial\phi_i}{\partial x}\phi_j
 *      +C_i D\frac{\partial\phi_i}{\partial x}\frac{\partial\phi_j}{\partial x}
 * \right\} dx
 * - \left[ D \frac{\partial C}{\partial x}\right]_0^L = 0
 * \f]
 */

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
#include "common.h"

extern choi_okos *comp_global;

/**
 * Function used to recalculate the Jacobian matrix at each time step. This
 * function is supplied to the FEM solver and integrated to give the J matrix.
 * It represents every term in the FEM equation being solved except for the
 * source term, the boundary conditions, and the time-dependent portion (those
 * terms containing \f$\frac{\partial T_i}{\partial t}\f$).
 * @param p Finite element structure
 * @param guess Estimated values for the solution at the current time step
 * @param elem Element to integrate over. Since nearly everything in this
 *      function is in terms of local coordinates, the specific element is
 *      required.
 * @param x Current local coordinate \f$\xi\f$ to calculate the residual at.
 *      This value is chosen by the integration function.
 * @param f1 Value for \f$i\f$ in \f$\phi_i\f$
 * @param f2 Value for \f$j\f$ in \f$\phi_j\f$
 * 
 * @returns Calculated value for the residual (not integrated)
 */
double ResVap(struct fe1d *p, matrix *guess, Elem1D *elem,
               double x, int f1, int f2)
{
           /* Used to store the values of the three terms of the PDE calculated
            * by this function */
    double term1 = 0, term2 = 0, term3 = 0,
           /* These hold the values for diffusivity and gradient of diffusivity
            * once they've been calculated */
           DDDx = 0, D = 0,
           /* D and grad(D) on the nodes */
           Di, Pi, Ti,
           /* Final value of the function */
           value = 0;
    int i;
    basis *b;
    b = p->b;

    /* Create a solution structure for the current guess so that the temperature
     * at any given point can be calculated. */
    solution *s;
    s = CreateSolution(p->t, p->dt, guess);

    /* Calculate thermal conductivity, density, heat capacity, and thermal
     * conductivity gradient at x */
    for(i=0; i<b->n; i++) {
        //Ti = EvalSoln1D(p, TVAR, elem, s, valV(elem->points, i));
        Pi = EvalSoln1D(p, CVAR, elem, s, valV(elem->points, i));
        //Di = DIFF(uscaleTemp(p->chardiff, Ci), uscaleTemp(p->charvals, Ti));
        Di = DIFF(uscaleTemp(p->chardiff, Pi), TINIT);
        printf("Di = %g\n", Di);
        D += Di * b->phi[i](x);
        DDDx += Di * b->dphi[i](x);
    }
    /* Then delete the temporary solution we made earlier. */
    free(s);

    /* Calculate the value of the term involving the second derivative of
     * temperature */
    term1 = D;
    term1 *= b->dphi[f1](x) * b->dphi[f2](x);
    term1 *= IMap1D(p, elem, x);

    /* Now that we have the gradient of thermal conductivity, we can calculate
     * the value of the term containing it. */
    term2 = DDDx * b->dphi[f1](x) * b->phi[f2](x);
    term2 *= 1/IMap1D(p, elem, x);
    /* Hopefully the above term is correct. It appears to yield accurate
     * results, at least. */

    /* This determines the value of the term that arises due to the moving
     * mesh. */
    term3 = b->dphi[f1](x) * IMapDt1D(p, elem, x);
    term3 *= b->phi[f2](x);
    term3 *= 1/IMap1D(p, elem, x);
    
    /* Combine all the terms and return the result */
    value = (term1 - term2)/p->chardiff.alpha + term3;
    return value;
}

/* Calculate the coefficient matrix for the time derivative unknowns. This does
 * the same thing as the "Residual" function, only it calculates the matrix
 * multiplying \f$\frac{\partial T_i}{\partial t}\f$.
 * @param p Finite element structure
 * @param guess Estimated values for the solution at the current time step
 * @param elem Element to integrate over. Since nearly everything in this
 *      function is in terms of local coordinates, the specific element is
 *      required.
 * @param x Current local coordinate \f$\xi\f$ to calculate the residual at.
 *      This value is chosen by the integration function.
 * @param f1 Value for \f$i\f$ in \f$\phi_i\f$
 * @param f2 Value for \f$j\f$ in \f$\phi_j\f$
 * 
 * @returns Calculated value for the residual (not integrated)
 * 
 * @see ResMass
 */
double ResDtVap(struct fe1d *p, matrix *guess, Elem1D *elem,
                 double x, int f1, int f2)
{
    double value;
    basis *b;
    b = p->b;

    value = b->phi[f1](x) * b->phi[f2](x) / IMap1D(p, elem, x);

    return value;
}

/**
 * Determine the external temperature used when heating.
 * @param p Finite element problem structure
 * @param row Node number to look at (Not used)
 * @returns Dimensionless external temperature
 */
double ExternalConcVap(struct fe1d *p, int row)
{
    return scaleTemp(p->chardiff, p->chardiff.Te);
    return 0;
}

/**
 * Calculate the value of \f$\underline{n}\cdot\nabla T\f$ on the boundary
 * @param p Finite element problem structure
 * @param row Row to calculate stuff at. This is guaranteed to be on the
 *      boundary by the solver.
 * @returns Temperature gradient at the surface
 */
double ConvBCVap(struct fe1d *p, int row)
{
    double Pinf = ExternalConc(p, row);
    double P;
    double Bi = BiotNumber(p->chardiff);

    /* Calculate the convective boundary condition based on the current estimate
     * for the temperature at the boundary. */
    if(p->guess) {
        P = val(p->guess, row, 0);
        return Bi*(P-Pinf);
    } else {
        return 0;
    }
}

