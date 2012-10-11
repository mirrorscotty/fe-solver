#ifndef MATRIX_H
#define MATRIX_H

/* Macro to shorten code where a number is added to the current value of an
 * element in a matrix. */
#define addval(A, VAL, I, J) setval((A), (VAL) + val((A), (I), (J)), (I), (J))

/**
 * @struct matrix
 * @brief A data structure for holding two-dimensional matricies
 * @var matrix::array
 * A pointer to the raw data
 * @var matrix::rows
 * Number of rows in the matrix
 * @var matrix::cols
 * Number of columns
 */
typedef struct {
    double **array;
    int rows;
    int cols;
} matrix;

void DestroyMatrix(matrix*);
matrix* CalcMinor(matrix*, int, int);
double CalcDeterminant(matrix*);
int nCols(matrix*);
int nRows(matrix*);
double val(matrix*, int, int);
void setval(matrix*, double, int, int);
matrix* CreateMatrix(int, int);
matrix* CreateOnesMatrix(int, int);
matrix* mtxtrn(matrix*);
matrix* mtxmul(matrix*, matrix*);
matrix* mtxmulconst(matrix*, double k);
matrix* mtxadd(matrix*, matrix*);
matrix* mtxneg(matrix*);
matrix* CalcAdj(matrix*);
matrix* CalcInv(matrix*);

matrix* ExtractColumn(matrix*, int);
matrix* AugmentMatrix(matrix*, matrix*);

matrix* ParseMatrix(char*);
matrix* linspace(double, double, int);
void Map(matrix*, double (*func)(double));
void mtxprnt(matrix*);
void mtxprntfile(matrix*, char*);

/**
 * @struct vector
 * @brief A data structure for a vector of abitrary length
 * @var vector::v
 * A pointer to the raw data
 * @var vector::length
 * The number of components in the vector
 */
typedef struct {
    double *v;
    int length;
} vector;

vector* CreateVector(int);
void DestroyVector(vector*);
vector* linspaceV(double, double, int);
vector* CopyVector(vector*);
//double valV(vector*, int);
void setvalV(vector*, int, double);
int len(vector*);
void PrintVector(vector*);

vector* addV(vector*, vector*);
vector* subtractV(vector*, vector*);
double dotV(vector*, vector*);
vector* scalarmultV(double, vector*);
int equalV(vector*, vector*);

#define valV(VECTOR, INDEX) (VECTOR)->v[(INDEX)]

#endif
