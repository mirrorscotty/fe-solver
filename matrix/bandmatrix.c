#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "bandmatrix.h"
#include "matrix.h"

bndmatrix* CreateBandMatrix(int rows, int bandwidth)
{
    int i;
    bndmatrix *bm;
    bm = (bndmatrix*) calloc(1, sizeof(bndmatrix));

    bm->r = rows;
    bm->w = bandwidth;

    bm->data = (double**) calloc(rows, sizeof(double*));
    for(i=0; i<rows; i++)
        bm->data[i] = (double*) calloc(bandwidth, sizeof(double));

    return bm;
}

void DestroyBandMatrix(bndmatrix *bm)
{
    int i;
    if(bm) {
        for(i=0; i<bm->r; i++)
            free(bm->data[i]);
        free(bm);
    }
}

double valB(bndmatrix *bm, int row, int col)
{
    int datacol;

    /* Determine if the indicies given are valid */
    if((row>=bm->r) || (col>=bm->r)) {
        fprintf(stderr, "Index out of bounds.\n");
        return 0;
    }

    /* Integer division truncates the result. */
    datacol = col+bm->w/2-row;

    /* Return 0 if the indicated element in the matrix is outside of the banded
     * region. */
    if((datacol<0) || (datacol>=bm->w))
        return 0;

    /* Otherwise, return the correct value. */
    return bm->data[row][datacol];
}

double setvalB(bndmatrix *bm, double val, int row, int col)
{
    int datacol;

    /* Determine if the indicies given are valid */
    if((row>=bm->r) || (col>=bm->r)) {
        fprintf(stderr, "Index out of bounds.\n");
        return 0;
    }

    /* Integer division truncates the result. */
    datacol = col+bm->w/2-row;

    /* Return 0 if the indicated element in the matrix is outside of the banded
     * region. */
    if((datacol<0) || (datacol>=bm->w)) {
        fprintf(stderr, "Cannot set value. Bandwidth too small.\n");
        return 0;
    }

    bm->data[row][datacol] = val;

    return val;
}

bndmatrix* ConvertFromDenseMatrix(matrix *source, int bw)
{
    int i, j;
    bndmatrix *dest;

    dest = CreateBandMatrix(nRows(source), bw);

    for(i=0; i<dest->r; i++) {
        for(j=0; j<bw; j++) {
            if((j+i-bw/2 >= 0) && (j+i-bw/2 < dest->r))
                setvalB(dest, val(source, i, j+i-bw/2), i, j+i-bw/2);
        }
    }

    return dest;
}

void bandprnt(bndmatrix *bm)
{
    int i, j;
    int datacol;
    double v;

    for(i=0; i<bm->r; i++) {
        printf("[ ");
        for(j=0; j<bm->r; j++) {
            datacol = j+bm->w/2-i;
            if((datacol<0) || (datacol>=bm->w)) {
                printf("%e ", 0.0);
            } else {
                v = valB(bm, i, j);
                if(fabs(v) < 1e-14)
                    printf("%e ", 0.0);
                else
                    printf("%e ", v);
            }
        }
        printf("]\n");
    }
}
