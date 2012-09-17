#ifndef BANDMATRIX_H
#define BANDMATRIX_H

#include "matrix.h"

typedef struct {
    double **data; /* Raw data */
    int r; /* Rows */
    int w; /* Bandwidth */
} bndmatrix;

bndmatrix* CreateBandMatrix(int, int);
void DestroyBandMatrix(bndmatrix*);
double valB(bndmatrix*, int, int);
double setvalB(bndmatrix*, double, int, int);
bndmatrix* ConvertFromDenseMatrix(matrix*, int);
void bandprnt(bndmatrix*);

#endif
