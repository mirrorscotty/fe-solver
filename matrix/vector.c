/**
 * @file vector.c
 * Define lots of functions for manipulating vectors of arbitrary length
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "matrix.h"

/**
 * @brief Create a 1d vector of length n
 * @param n Number of elements in the vector
 * @return The newly created vector
 */
vector* CreateVector(int n)
{
    vector *v;
    v = (vector*) calloc(1, sizeof(vector));

    v->v = (double*) calloc(n, sizeof(double));
    v->length = n;

    return v;
}

/**
 * @brief Free the memory allocated for the vector.
 * @param The pointer to the vector to deallocate memory for
 */
void DestroyVector(vector *v)
{
    free(v->v);
    free(v);
    return;
}

/**
 * @brief Print the vector to stdout
 * @param A pointer to the vector to print
 */
void PrintVector(vector *v)
{
    int i;

    printf("[");
    for(i=0; i<len(v)-1; i++) {
        printf(" %g,", valV(v, i));
    }
    printf(" %g ]\n", valV(v, len(v)-1));
}

/**
 * @brief Make a vector of equally spaced points.
 * @param start The value of the first comonent of the vector
 * @param end The value of the last
 * @param nelem The total number of comonents
 * @returns The vector of points
 *
 * @see linspace
 */
vector* linspaceV(double start, double end, int nelem)
{
    vector *x;
    int i;

    x = CreateVector(nelem);
    for(i=0; i<nelem; i++) {
        setvalV(x, i, start + i*(end-start)/(nelem-1));
    }

    return x;
}

/**
 * @brief Copy a vector
 * @param v A pointer to the original vector
 * @returns A pointer to a new vector with the same components as v
 */
vector* CopyVector(vector *v)
{
    vector *new;
    int i;
    new = CreateVector(len(v));
    for(i=0; i<len(v); i++) {
        setvalV(new, i, valV(v, i));
    }
    return new;
}

/* Return the i-th value */
//double valV(vector *v, int i)
//{
//    if(i<0 || i>=v->length) {
//        /* Index out of bounds */
//        return 0;
//    }
//
//    return v->v[i];
//}
/* Replaced by a macro for efficiency */

/**
 * @brief Set the i-th value to "val"
 * @param v The vector to set the value in
 * @param i The index of the component to change
 * @param val The new value the component should have
 */
void setvalV(vector *v, int i, double val)
{
    if(i<0 || i>=v->length) {
        /* Index out of bounds */
        return;
    }

    v->v[i] = val;
    return;
}

/**
 * @brief Return the length of the vector
 * @param v The vector to get the length for
 * @returns The length of the vector
 */
int len(vector* v)
{
    return v->length;
}

/**
 * Add two vectors together, element by element
 *
 * The two vectors should be of equal length
 *
 * @param a The first vector
 * @param b Second vector
 * @returns a+b
 */
vector* addV(vector *a, vector *b)
{
    int i;
    vector *c;
    c = CreateVector(len(a));

    for(i=0; i<len(a); i++) {
        setvalV(c, i, valV(a, i) + valV(b, i));
    }
    return c;
}

/**
 * Subtract vector b from vector a 
 *
 * @param a The first vector
 * @param b The vector subtracted from a
 * @returns a-b
 */
vector* subtractV(vector *a, vector *b)
{
    int i;
    vector *c;
    c = CreateVector(len(a));

    for(i=0; i<len(a); i++) {
        setvalV(c, i, valV(a, i) - valV(b, i));
    }
    return c;
}

/**
 * @brief Calculate the dot product of two vectors
 * @param a A vector of arbitrary length
 * @param b A vector of the same length
 * @returns a dot b
 */
double dotV(vector *a, vector *b)
{
    int i;
    double result = 0;

    for(i=0; i<len(a); i++) {
        result += valV(a, i) * valV(b, i);
    }
    return result;
}

/**
 * @brief Multiply a vector by a scalar
 * @param k The constant to multiply each component by
 * @param v The vector to multiply
 * @returns k*v
 */
vector* scalarmultV(double k, vector *v)
{
    int i;
    vector *c;

    c = CreateVector(v->length);

    for(i=0; i<len(v); i++) {
        setvalV(c, i, k*valV(v, i));
    }

    return c;
}

/**
 * @brief Determine if vectors are equal
 * 
 * Compare two vectors element by element to determine equality
 *
 * @param a The first vector
 * @param b Second vector
 * @returns 1 if equal, 0 if not equal
 */
int equalV(vector *a, vector *b)
{
    int i;
    double tol = 1e-10;

    if(a == b)
        return 1;

    if(len(a) != len(b)) {
        return 0;
    }

    for(i=0; i<len(a); i++) {
        if(fabs(valV(a, i) - valV(b, i)) >= tol)
            return 0;
    }

    return 1;
}

