#include <stdio.h>
#include <stdlib.h>

#include "matrix.h"

/* Create a 1d vector of length n */
vector* CreateVector(int n)
{
    vector *v;
    v = (vector*) calloc(1, sizeof(vector));

    v->v = (double*) calloc(n, sizeof(double));
    v->length = n;

    return v;
}

/* Free the memory allocated for the vector. */
void DestroyVector(vector *v)
{
    free(v->v);
    free(v);
    return;
}

/* Print the vector to stdout */
void PrintVector(vector *v)
{
    int i;

    printf("[");
    for(i=0; i<len(v)-1; i++) {
        printf(" %g,", valV(v, i));
    }
    printf(" %g ]\n", valV(v, len(v)-1));
}

/* Return the i-th value */
double valV(vector *v, int i)
{
    if(i<0 || i>=v->length) {
        /* Index out of bounds */
        return 0;
    }

    return v->v[i];
}

/* Set the i-th value to "val" */
void setvalV(vector *v, int i, double val)
{
    if(i<0 || i>=v->length) {
        /* Index out of bounds */
        return;
    }

    v->v[i] = val;
    return;
}

/* Return the length of the vector */
int len(vector* v)
{
    return v->length;
}

/* Add two vectors together */
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

/* Subtract vector b from vector a */
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

/* Return the dot product of two vectors */
double dotV(vector *a, vector *b)
{
    int i;
    double result = 0;

    for(i=0; i<len(a); i++) {
        result += valV(a, i) * valV(b, i);
    }
    return result;
}

/* Multiply a vector by a scalar */
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

