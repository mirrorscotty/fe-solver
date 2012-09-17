#include "matrix.h"
#include "mtxsolver.h"

void SwapRows(matrix*, int, int);
int FindPivot(matrix*, int);

/* Code from http://compprog.wordpress.com/2007/12/11/gaussian-elimination/ */

void ForwardSubstitution(matrix* a) {
    int i, j, k, max;
    int n = nRows(a);
    double t;
    for (i = 0; i < n; ++i) {
        max = i;
        for (j = i + 1; j < n; ++j)
            if (fabs(val(a, j, i)) > fabs(val(a, max, i)))
                max = j;

        for (j = 0; j < n + 1; ++j) {
            t = val(a, max, j);
            setval(a, val(a, i, j), max, j);
            setval(a, t, i, j);
        }

        for (j = n; j >= i; --j)
            for (k = i + 1; k < n; ++k)
                setval(a, val(a, k, j) - val(a, k, i)/val(a, i, i) * val(a, i, j), k, j);

    }
}

void ReverseElimination(matrix *a) {
    int i, j;
    int n = nRows(a);
    for (i = n - 1; i >= 0; --i) {
        setval(a, val(a, i, n)/val(a, i, i), i, n);
        setval(a, 1, i, i);
        for (j = i - 1; j >= 0; --j) {
            setval(a, val(a, j, n) - val(a, j, i) * val(a, i, n), j, n);
            setval(a, 0, j, i);
        }
    }
}

/* Solve an matrix equation of the form Ax=B, where x is the vector of
 * unknowns. */
matrix* SolveMatrixEquation(matrix *A, matrix *B)
{
    matrix *C, *u;
    C = AugmentMatrix(A, B);
    ForwardSubstitution(C);
    ReverseElimination(C);
    u = ExtractColumn(C, nCols(A));
    DestroyMatrix(C);
    return u;
}
