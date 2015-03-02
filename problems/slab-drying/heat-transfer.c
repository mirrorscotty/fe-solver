/**
 * @file heat-transfer.c
 * Solves the nonlinear heat conduction equation. 
 * \f[
 * \rho C_p\frac{\partial T}{\partial t} = \nabla\cdot\left(k\nabla T\right)
 * \f]
 * Density, heat capacity, and thermal conductivity are all temperature-
 * dependent and are calculated using the Choi-Okos equations. Because density
 * is changing, but mass remains constant, there is a certain degree of
 * shrinkage/expansion that occurs during heating. This is calculated from the
 * deformation gradient as follows:
 * \f[
 * \underline{\underline{F}}(\underline{X}, t) = \frac{\rho(0)}{\rho(t)}
 * \f]
 * The updated nodal values are calculated at each time step using this tensor
 * and the mesh velocity is incorporated into the FEM equations.
 *
 * Finite Element Formulation (Global Coordinates):
 * \f[
 * R^i = \int_0^L \left\{
 *      \rho C_p \frac{\partial T_i}{\partial t}\phi_i\phi_j
 *      + T_i \frac{\partial\phi_i}{\partial t}\phi_j
 *      -T_i\frac{\partial k}{\partial x}\frac{\partial\phi_i}{\partial x}\phi_j
 *      +T_i k\frac{\partial\phi_i}{\partial x}\frac{\partial\phi_j}{\partial x}
 * \right\} dx
 * - \left[ k \frac{\partial T}{\partial x}\right]_0^L = 0
 * \f]
 * Here, the first term is the rate of change of temperature at a single point
 * in the slab. The second term relates to the rate of deformation of the solid
 * and how that affects it's temperature. The third term accounts for the
 * extra derivative required due to the temperature dependence of thermal
 * diffusivity, and the final term in the integral is the \f$\nabla^2 T\f$ term.
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

#define CONST_MC 0

/* The following two functions are for heat conduction. */
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
double ResHeat(struct fe1d *p, matrix *guess, Elem1D *elem,
		double x, int f1, int f2)
{
           /* Used to store the values of the three terms of the PDE calculated
            * by this function */
    double term1 = 0, term2 = 0, term3 = 0,
           /* These hold the values of the thermal properties once they've been
            * calculated */
           DkDx = 0, kval = 0, Cpval = 0, rhoval = 0,
           /* Values of the thermal properties on the nodes */
           Cpi, rhoi, Ti, Ci = CONST_MC, ki,
           /* Final value of the function */
           value = 0;
    int i;
    basis *b;
    b = p->b;
    choi_okos *coweti;

    /* Create a solution structure for the current guess so that the temperature
     * at any given point can be calculated. */
    solution *s;
    s = CreateSolution(p->t, p->dt, guess);

    /* Calculate thermal conductivity, density, heat capacity, and thermal
     * conductivity gradient at x */
    for(i=0; i<b->n; i++) {
        Ti = EvalSoln1D(p, TVAR, elem, s, valV(elem->points, i));
        //Ci = EvalSoln1D(p, CVAR, elem, s, valV(elem->points, i));
        //coweti = AddDryBasis(comp_global, Ci);
        coweti = comp_global;
        
        ki = k(coweti, uscaleTemp(p->charvals, Ti));
        kval += ki * b->phi[i](x);
        DkDx += ki * b->dphi[i](x);
        rhoi = rho(coweti, uscaleTemp(p->charvals, Ti));
        rhoval += rhoi * b->phi[i](x);
        Cpi = Cp(coweti, uscaleTemp(p->charvals, Ti));
        Cpval += Cpi * b->phi[i](x);

        //DestroyChoiOkos(coweti);
    }
    /* Then delete the temporary solution we made earlier. */
    free(s);

    /* Calculate the value of the term involving the second derivative of
     * temperature */
    term1 = kval;
    term1 *= b->dphi[f1](x) * b->dphi[f2](x);
    term1 *= IMap1D(p, elem, x);

    /* Now that we have the gradient of thermal conductivity, we can calculate
     * the value of the term containing it. */
    term2 = DkDx * b->dphi[f1](x) * b->phi[f2](x);
    term2 *= 1/IMap1D(p, elem, x);

    /* This determines the value of the term that arises due to the moving
     * mesh. */
    term3 = b->dphi[f1](x) * IMapDt1D(p, elem, x);
    term3 *= b->phi[f2](x);
    term3 *= 1/IMap1D(p, elem, x);
    
    /* Combine all the terms and return the result */
    value = (term1 - term2)/(rhoval*Cpval)/p->charvals.alpha + term3;
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
 * @see ResHeat
 */
double ResDtHeat(struct fe1d *p, matrix *guess, Elem1D *elem,
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
double ExternalTemp(struct fe1d *p, int row)
{
    return scaleTemp(p->charvals, p->charvals.Te);
}

/**
 * Calculate the value of \f$\underline{n}\cdot\nabla T\f$ on the boundary
 * @param p Finite element problem structure
 * @param row Row to calculate stuff at. This is guaranteed to be on the
 *      boundary by the solver.
 * @returns Temperature gradient at the surface
 */
double ConvBCHeat(struct fe1d *p, int row)
{
    double Tinf = ExternalTemp(p, row);
    double T;
    double Bi = BiotNumber(p->charvals);

    /* Calculate the convective boundary condition based on the current estimate
     * for the temperature at the boundary. */
    if(p->guess) {
        T = FetchGuessValue(p, row, TVAR);
        printf("T = %g, Tinf = %g\n", T, Tinf);
        return -Bi*(T-Tinf);
    } else {
        return 0;
    }
}

