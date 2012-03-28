#include "regress.h"
#include "matrix.h"

/* Equivalent of the Matlab "regress" function. */
matrix* regress(matrix *y, matrix *X)
{
    matrix *Xt, *XtX, *XtXinv, *XtXinvXt, *beta;

    Xt = mtxtrn(X);
    XtX = mtxmul(Xt, X);
    XtXinv = CalcInv(XtX);
    XtXinvXt = mtxmul(XtXinv, Xt);

    beta = mtxmul(XtXinvXt, y);

    DestroyMatrix(Xt);
    DestroyMatrix(XtX);
    DestroyMatrix(XtXinv);
    DestroyMatrix(XtXinvXt);

    return beta;
}

/* Matlab "polyfit" function. */
matrix* polyfit(matrix* x, matrix* y, int order)
{
    matrix *Y, *X, *beta;
    int i, j, nelem;

    X = NULL;
    nelem = mtxlen1(x); /* Row matrix */

    X = CreateMatrix(nelem, order+1);

    for(i=0; i<=order; i++) {
        for(j=0; j<nelem; j++) {
            setval(X, pow(val(x, 0, j), i), j, i);
        }
    }
    Y = mtxtrn(y);

    beta = regress(Y, X);
    DestroyMatrix(X);
    DestroyMatrix(Y);

    return beta;
}
 
