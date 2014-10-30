/**
 * @file isoparam.c
 * All sorts of stuff for isoparametric mapping!
 */

#include <stdio.h>
#include "isoparam.h"

/**
 * @brief Macro to simplify the function definitions below.
 *
 * Otherwise, it'd be nearly
 * the same thing copy and pasted 4 times. No one wants that!
 * These functions handle isoparametric mapping for both bilinear and
 * biquad elements.
 * Currently, the Elem2D structure can't handle anything but bilinear elements,
 * so this won't work with anything else either.
 * Also, this whole idea is probably a bad one.
 */
#define IMAP(NAME, GVAR, LVAR) \
    double IMap##NAME(struct fe *p, Elem2D *elem, double xi, double eta) \
    { \
        basis *b; \
        int i; \
        b = p->b; \
        int n = b->n; \
        double result = 0; \
\
\
        for(i=0; i<n; i++) { \
            result += b->Eval2D##GVAR(b, i, xi, eta) \
                * valV(elem->points[i], LVAR); \
        } \
\
        return result; \
    }

/* dx/dxi */
IMAP(XXi, x, 0)
/* dy/dxi */
IMAP(YXi, x, 1)
/* dx/deta */
IMAP(XEta, y, 0)
/* dy/deta */
IMAP(YEta, y, 1)

/**
 * @brief Calculate the Jacobian
 * @param p The problem definition
 * @param elem The element of interest
 * @param xi The local x-coordinate to evaluate at
 * @param eta The local y-coordinate
 * @returns Value of the Jacobian
 */
double IMapJ(struct fe *p, Elem2D *elem, double xi, double eta)
{
    return IMapXXi(p, elem, xi, eta)
        * IMapYEta(p, elem, xi, eta)
        - IMapYXi(p, elem, xi, eta)
        * IMapXEta(p, elem, xi, eta);
}

/**
 * @brief Calculate the differential volume for use when integrating in
 * cylindrical coordinates
 * @param p Problem structure
 * @param elem Element where stuff is being calculated
 * @param xi Local x-coordinate
 * @param eta Local y-coordinate
 * @returns Differential volume at that element
 */
double IMapCyl(struct fe *p, Elem2D *elem, double xi, double eta)
{
    double x = 0;
    basis *b;
    int n, i;

    b = p->b;
    n = b->n;

    for(i=0; i<n; i++)
        x += b->Eval2D(b, i, xi, eta) * valV(elem->points[i], 0);

    //return 1;
    return x;
}

/**
 * Same thing as the 2D version, only for one dimension.
 * @see IMapCyl
 */
double IMapCyl1D(struct fe1d *p, Elem1D *elem, double xi)
{
    double x = 0;
    basis *b;
    int n, i;

    b = p->b;
    n = b->n;

    for(i=0; i<n; i++)
        x += b->phi[i](xi) * valV(elem->points, i);

    return x;
}

double IEvalLin2Dx(struct fe *p, Elem2D *elem, int func, double x, double y)
{
    double r;
    basis *b;
    b = p->b;

    r = 0;
    r = EvalLin2Dx(b, func, x, y)*IMapYEta(p, elem, x, y);
    r -= EvalLin2Dy(b, func, x, y)*IMapYXi(p, elem, x, y);
    r *= IMapJ(p, elem, x, y);
    return r;
}

double IEvalLin2Dy(struct fe *p, Elem2D *elem, int func, double x, double y)
{
    double r;
    basis *b;
    b = p->b;

    r = 0;
    r -= EvalLin2Dx(b, func, x, y)*IMapXEta(p, elem, x, y);
    r += EvalLin2Dy(b, func, x, y)*IMapXXi(p, elem, x, y);
    r *= IMapJ(p, elem, x, y);
    return r;
}

/**
 * @brief Do isoparametric mapping for 1D problems.
 *
 * This function returns the value of the derivative of the global x-coordinate
 * with respect to the local x-coordinate.
 * \f[\frac{\partial x}{\partial \xi}\f]
 * @param p The problem being solved
 * @param elem The element to calculate stuff for
 * @param xi Local x-coordinate
 * */
double IMap1D(struct fe1d *p, Elem1D *elem, double xi)
{
    /* Only works for linear elements */
    double x1 = valV(elem->points, 0);
    double x2 = valV(elem->points, 1);
    
    double result = 0;
    basis *b;
    b = p->b;
    int n = p->b->n;
    int i;

    for(i=0; i<n; i++) {
        result += b->phi[i](xi) * valV(elem->points, i);
    }

    /* Cheat and return the right value. This only works reliably for linear
     * interpolation functions */
    result = 1/(x2-x1);

    return result;
}

/**
 * Calculate mesh velocity. (At least that's kinda what this is...)
 * \f[
 * v = \frac{\partial}{\partial t}\left[\frac{\partial \xi}{\partial x}\right]
 * \f]
 * @param p Finite element problem structure
 * @param elem Element to operate on
 * @param xi Local x-coordinate
 * @returns Value of the derivative above.
 */
double IMapDt1D(struct fe1d *p, Elem1D *elem, double xi)
{
    double DxiDx0, DxiDx1;

    if(!(elem->prev))
        return 0;

    DxiDx0 = 1/IMap1D(p, elem, xi); /* DxiDx(t=n) */
    DxiDx1 = 1/IMap1D(p, elem->prev, xi); /* DxiDx(t=n-1) */

    /* Just dividing by dt is a lazy way to do this. The correct way would be to
     * find the time step size for the current step and divide by that. Since
     * the solver currently doesn't change time steps, this isn't important, and
     * the formula here works fine. */
    return (DxiDx0-DxiDx1)/p->dt;
}

