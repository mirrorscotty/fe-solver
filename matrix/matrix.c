/**
 * @file matrix.c
 * This file contains all sorts of nifty matrix functions. All of the matricies
 * here have indicies that start with 0 instead of 1.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "matrix.h"

/**
 * @brief Create a matrix of equally spaced points
 *
 * This is the equivalent of the Matlab linspace command.
 * 
 * @param start First value
 * @param end Last value
 * @param Number of points to put in the matrix
 * @return Row matrix containing "nelem" points
 */
matrix* linspace(double start, double end, int nelem)
{
    matrix* x;
    int i;
    
    x = CreateMatrix(1, nelem);
    for(i=0; i<nelem; i++) {
        setval(x, start + i*(end-start)/(nelem-1), 0, i);
    }

    return x;
}
    
/**
 * @brief Get the number of columns in a 1D matrix of doubles
 * @param A The matrix!
 * @return Number of columns 
 */
int nCols(matrix *A)
{
    return A->cols;
}

/**
 * @brief Return the number of rows in a two dimensional matrix of doubles
 *
 * @param A The matrix of interest
 * @return Number of rows
 */
int nRows(matrix *A)
{
    return A->rows;
}

/**
 * @brief Get the value of a specific element in a matrix
 *
 * @param A The matrix containing the values
 * @param row The row to pull the value from
 * @param col The column to pull the value from
 * @return The value stored at row, col in matrix A
 */
double val(matrix *A, int row, int col)
{
    if((row >= nRows(A)) ||  (col >= nCols(A)) || (row < 0) || (col < 0)) {
        fprintf(stderr, "Error: index out of bounds. (%d, %d)\n", row, col);
        return NAN;
    }
    return A->array[row][col];
}

/**
 * Set the value of a particular element of a matrix
 *
 * @param A The matrix to set the value in
 * @param value The value to set
 * @param row The row
 * @param col Column
 */
void setval(matrix *A, double value, int row, int col)
{
    if((row >= nRows(A)) ||  (col >= nCols(A)) || (row < 0) || (col < 0)) {
        //fprintf(stderr, "Error: index out of bounds. (%d, %d)\n", row, col);
        return;
    }
    A->array[row][col] = value;
}

/**
 * @brief Print out a matrix
 * @param A The matrix to print
 */
void mtxprnt(matrix *A)
{
    int i, j;
    double v;
    
    for(i=0; i<nRows(A); i++) {
        printf("[ ");
        for(j=0; j<nCols(A); j++) {
            v = val(A, i, j);
            /* If the value is annoyingly close to zero, just print out zero
             * instead. */
            if(fabs(v) < 1e-14)
                printf("%e ", 0.0);
            else
                printf("%e ", v);
        }
        printf("]\n");
    }
}

/**
 * @brief Print a matrix out to some random file
 * @param A The matrix to print
 * @param filename The filename to print to
 */
void mtxprntfile(matrix *A, char *filename)
{
    int i, j;
    FILE *file;

    file = fopen(filename, "w");
    
    for(i=0; i<nRows(A); i++) {
        for(j=0; j<nCols(A); j++) {
            fprintf(file, "%e,", val(A, i, j));
        }
        fprintf(file, "\n");
    }
    fclose(file);
}

/**
 * @brief Transpose a matrix
 * @param x The matrix to transpose
 * @return The transpose of matrix x
 */
matrix* mtxtrn(matrix *x)
{
    matrix *xt;
    int rows = nRows(x);
    int cols = nCols(x);
    int i, j;

    xt = NULL;
/*
    if(nRows(x) != nCols(x)) {
        fprintf(stderr, "ERROR!");
        return xt;
    }
*/

    xt = CreateMatrix(cols, rows);

    for(i=0; i<cols; i++) {
        for(j=0; j<rows; j++) {
            setval(xt, val(x, j, i), i, j);
        }
    }

    return xt;
}

/**
 * @brief Multiply matricies using nifty index notation!
 *
 * Simply calculates A*B
 *
 * @param A The first matrix to multiply
 * @param B The second one!
 * @return A*B
 */
matrix* mtxmul(matrix *A, matrix *B)
{
    int Ar, Ac, Br, Bc;
    int i, j, k;
    matrix *C;

    C = NULL;
    
    Ar = nRows(A);
    Ac = nCols(A);
    Br = nRows(B);
    Bc = nCols(B);

    /* If the matricies dimensions aren't correct, return NULL */
    if(Ac != Br) {
        fprintf(stderr, "Error: Incompatible matrix dimensions.\n");
        return C;
    }
    
    /* Allocate Memory */
    C = CreateMatrix(Ar, Bc);
    
    /* Cik = AijBjk */
    for(i=0; i<Ar; i++) {
        for(k=0; k<Bc; k++) {
            for(j=0; j<Ac; j++) {
                C->array[i][k] += A->array[i][j] * B->array[j][k];
            }
        }
    }

    return C;
}

/**
 * @brief Multiply a matrix by a constant
 * This multiplies each element of a matrix by a constant.
 * @param A The matrix to multiply
 * @param k The scalar
 * @return The new matrix
 */
matrix* mtxmulconst(matrix *A, double k)
{
    int Ar, Ac;
    int i, j;
    matrix *C;

    Ar = nRows(A);
    Ac = nCols(A);

    C = CreateMatrix(Ar, Ac);

    for(i=0; i<Ar; i++) {
        for(j=0; j<Ac; j++) {
            setval(C, k*val(A, i, j), i, j);
        }
    }

    return C;
}

/**
 * @brief Add two matricies.
 * Todo: Add error checking to make sure the dimensions agree.
 * @param A Some random matrix
 * @param B Another random matrix with the same dimensions as A
 * @return A+B
 */
matrix* mtxadd(matrix *A, matrix *B)
{
    int rows = nRows(A);
    int cols = nCols(A);
    
    matrix *C;
    
    int i, j;
    
    C = CreateMatrix(rows, cols);
    
    for(i=0; i<rows; i++) {
        for(j=0; j<cols; j++) {
            setval(C, val(A, i, j) + val(B, i, j), i, j);
        }
    }
    
    return C;
}

/**
 * @brief Multiply a matrix by -1 in place.
 * @param A The matrix to multiply by -1
 * @return A pointer to the same matrix A that was passed to this function*/
matrix* mtxneg(matrix *A)
{
    int i, j;
    for(i=0; i<nRows(A); i++) {
        for(j=0; j<nCols(A); j++) {
            setval(A, -1*val(A, i, j), i, j);
        }
    }
    return A;
}

/**
 * @brief Calculate the minor of the row,col element of a matrix.
 *
 * Code here courtesy of http://www.daniweb.com/software-development/c/code/216687
 * @param A The matrix to calculate stuff for
 * @param row The row of the element
 * @param col The column the element is in
 * @return The value of the minor
 */
matrix* CalcMinor(matrix* A, int row, int col) {
    int i, j, a, b;
    int order;
    matrix *minor;

    minor = NULL;
    order = nRows(A);

    if(order <= 1)
        return NULL;
    if(row >= order || col >= order)
        return NULL;
    if( !(minor = CreateMatrix(order-1, order-1)) )
        return NULL;

    a = b = 0;

    for(i=0; i<order; i++) {
        if(i != row) {
            b = 0;
            for(j=0; j<order; j++) {
                if(j != col) {
                    setval(minor, val(A, i, j), a, b);
                    b++;
                }
            }
        a++;
        }
    }

    return minor;
}

/**
 * @brief Calculate the determinant of a matrix
 *
 * Borrowed from the same site as CalcMinor
 * @param p The matrix to calculate the determiant of. Must be square.
 * @return The determinant of p
 * @see CalcMinor
 */
double CalcDeterminant(matrix *p)
{
    int i, order;
    double result;
    matrix *minor;

    result = 0;
    minor = NULL;
    order = nRows(p);

    if(order < 1) {
        fprintf(stderr, "CalcDeterminant(): Invalid Matrix.");
        return 0;
    }

    if(order == 1)
        return val(p, 0, 0);

    for(i=0; i<order; i++) {
        if( !(minor = CalcMinor(p, 0, i)) ) {
            fprintf(stderr, "CalcDeterminant(): Memory allocation failed.");
            return 0;
        }

        result += ( pow(-1, i) * p->array[0][i] * CalcDeterminant(minor));

        DestroyMatrix(minor);
    }

    return result;
}

/**
 * @brief Calculate the adjugate matrix of A
 * 
 * @param A The matrix of interest
 * @return The adjugate matrix
 */
matrix* CalcAdj(matrix* A)
{
    int i, j;
    double cofactor;
    matrix *minor, *adj, *adjt;
    
    adj = CreateMatrix(nRows(A), nRows(A));
    for(i=0; i<nRows(A); i++) {
        for(j=0; j<nRows(A); j++) {
            minor = CalcMinor(A, i, j);
            cofactor = pow(-1, (i+j+2)) * CalcDeterminant(minor);
            DestroyMatrix(minor);
            setval(adj, cofactor, i, j);
        }
    }

    adjt = mtxtrn(adj);
    DestroyMatrix(adj);
    

    return adjt;
}

/**
 * @brief Calculate the inverse of A
 * This is a pretty slow algorithm and only works well for small matricies
 * @param A The original matrix
 * @returns The inverse of A
 */
matrix* CalcInv(matrix* A)
{
    matrix *inv, *adj;
    int i, j;
    double det;

    inv = NULL;
    det = CalcDeterminant(A);
    adj = CalcAdj(A);
    inv = CreateMatrix(nRows(A), nRows(A));
    
    for(i=0; i<nRows(A); i++) {
        for(j=0; j<nRows(A); j++) {
            setval(inv, val(adj, i, j)/det, i, j);
        }
    }

    return inv;
}

/**
 * @brief Make an agumented matrix from the two arguments
 * This is used to set up a matrix for the equation solver. Basically, matrix
 * B is stuck on the end of A like this: [A | B].
 * @param A The first matrix
 * @param B Second matrix
 * @returns A matrix composed of the first matrix on the left and the second
 * on the right.
 */
matrix* AugmentMatrix(matrix *A, matrix *B)
{
    int i, j;
    matrix *C;
    C = CreateMatrix(nRows(A), nCols(A)+nCols(B));
    for(i=0; i<nRows(A); i++) {
        for(j=0; j<nCols(A); j++) {
            setval(C, val(A, i, j), i, j);
        }
        for(j=0; j<nCols(B); j++) {
            setval(C, val(B, i, j), i, j+nCols(A));
        }
    }

    return C;
}

/**
 * @brief Pull out a column of a matrix.
 * @param A The matrix to pull the column from
 * @param col The column to get
 * @returns A column matrix containing the values from the desired column of A
 */
matrix* ExtractColumn(matrix* A, int col)
{
    matrix *B;
    int i;
    B = CreateMatrix(nRows(A), 1);

    for(i=0; i<nRows(A); i++) {
	setval(B, val(A, i, col), i, 0);
    }

    return B;
}

/**
 * @brief Apply a function to every value in a matrix
 *
 * No value is returned because the original matrix is operated on.
 *
 * @param A The matrix of values
 * @param func The function to apply. It should accept a single double as an
 * argument and return a double.
 */
void Map(matrix* A, double (*func)(double))
{
    int i, j;
    for(i=0; i<nRows(A); i++) {
        for(j=0; j<nCols(A); j++) {
            setval(A, (*func)(val(A, i, j)), i, j);
        }
    }
}

/**
 * @brief Make a matrix!
 *
 * Create a matrix of the specified dimensions and initialize all the values
 * and initialize all the values to zero.
 * 
 * @param row The number of rows
 * @param col Number of columns
 * @returns The new matrix
 */
matrix* CreateMatrix(int row, int col)
{
    matrix *A;
    int i;
    A = NULL;

    if((row == 0) || (col == 0)) {
        printf("Matrix too small.");
        return A;
    }

    A = (matrix*) malloc(sizeof(matrix));
    if(!A) {
        fprintf(stderr, "Memory allocation error: %s\n", strerror(errno));
    }

    A->array = NULL;
    A->rows = 0;
    A->cols = 0;

    A->array = (double**) calloc(row, sizeof(double*));
    if(!A->array) {
        fprintf(stderr, "Memory allocation error: %s\n", strerror(errno));
        return A;
    }

    for(i=0; i<row; i++) {
        A->array[i] = (double*) calloc(col, sizeof(double));
        if(!A->array[i]) {
            fprintf(stderr, "Failed to allocate %lu bytes: %s\n", (unsigned long) sizeof(double)*col, strerror(errno));
            return A;
        }
    }

    A->rows = row;
    A->cols = col;
    
    return A;
}

/**
 * @brief Make a matrix of ones
 * 
 * Basically does the same thing as CreateMatrix, only it initializes each value
 * to one instead of zero.
 *
 * @param row Number of rows
 * @param col Number of columns
 * @returns A matrix
 */
matrix* CreateOnesMatrix(int rows, int cols)
{
    int i, j;
    matrix *A;
    A = CreateMatrix(rows, cols);
    for(i=0; i<rows; i++) {
        for(j=0; j<cols; j++) {
            setval(A, 1, i, j);
        }
    }
    return A;
}

///Maximum number of rows to allow
#define MAXROWS 100
///Maximum number of columns
#define MAXCOLS 100
///The length of each line of text that is parsed
#define LINELENGTH 80
/**
 * @brief Create a matrix from a line of text.
 *
 * This takes a random string of characters of the format [4,5,2;5,1,4] and
 * turns it into a nifty matrix.
 *
 * @param raw The character pointer to the string of characters
 * @returns A matrix created from the input
 */
matrix* ParseMatrix(char* raw)
{
    char *processed, *tmp;
    char **rows;
    double **values;
    matrix *out;
    int i, j;
    int nrows, ncols;

    processed = calloc(sizeof(char), MAXROWS*LINELENGTH);
    tmp = processed;

    while(*raw) {
        if(*raw == '[')
            *tmp = ' ';
        else if(*raw == ']')
            *tmp = ' ';
        else
            *tmp = *raw;
        tmp++;
        raw++;
    }

    /* Allocate memory to store the rows of the matrix as strings */
    values = (double**) calloc(sizeof(double*), MAXROWS);
    rows = (char**) malloc(sizeof(char*)*MAXROWS);
    for(i=0; i<MAXROWS; i++) {
        values[i] = (double*) calloc(sizeof(double), MAXCOLS);
        rows[i] = (char*) calloc(sizeof(char), LINELENGTH);
    }

    /* Row 1: */
    i = 0;
    tmp = strtok(processed, ";");
    strncpy(rows[i], tmp, LINELENGTH);

    /* All subsequent rows */
    for(i=1; (tmp = strtok(NULL, ";")); i++) {
        strncpy(rows[i], tmp, LINELENGTH);
    }

    nrows = i;

    for(i=0; i<nrows; i++) {
        /* First column */
        j=0;
        tmp = strtok(rows[i], ",");
        values[i][j] = atof(tmp);

        ncols = 1;
        /* Rest of the columns */
        for(j=1; (tmp = strtok(NULL, ",")); j++) {
            if(j>ncols-1)
                ncols = j+1;
            values[i][j] = atof(tmp);
        }
    }

    out = CreateMatrix(nrows, ncols);
    for(i=0; i<nrows; i++) {
        for(j=0; j<ncols; j++)
            setval(out, values[i][j], i, j);
    }

    /* Free stuff */
    free(processed);
    for(i=0; i<MAXROWS; i++) {
        free(values[i]);
        free(rows[i]);
    }
    free(rows);
    free(values);
    
    return out;
}


/**
 * @brief Free the memory allocated by CreateMatrix
 * @param A The matrix to destroy
 */
void DestroyMatrix(matrix *A)
{
    int i;

    for(i=0; i<A->rows; i++) {
        free(A->array[i]);
    }
    free(A->array);
    free(A);
}

