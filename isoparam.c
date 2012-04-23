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
        int n = b->n; \
        double result = 0; \
\
        b = p->b; \
\
        for(i=0; i<n; i++) { \
            result += EvalLin2D##GVAR(b, i, xi, eta) \
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

