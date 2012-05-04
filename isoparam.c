#include <stdio.h>
#include "isoparam.h"

/* Macro to simplify the function definitions below. Otherwise, it'd be nearly
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

/* Calculate the Jacobian */
double IMapJ(struct fe *p, Elem2D *elem, double xi, double eta)
{
    return IMapXXi(p, elem, xi, eta)
        * IMapYEta(p, elem, xi, eta)
        - IMapYXi(p, elem, xi, eta)
        * IMapXEta(p, elem, xi, eta);
}

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

