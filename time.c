#include "time.h"
#include "solution.h"
#include "matrix.h"

double DTime1D(struct fe1d* p, Elem1D *elem, double xi)
{
    solution *cur;
    double result = 0;
    int t = p->t-1; /* Time index for the most recent solution */

    cur = FetchSolution(p, t);
    //prev = FetchSolution(p, t-1);

    result = IEvalSoln1D(p, elem, cur, xi);// - IEvalSoln1D(p, elem, prev, xi);
    //result = result/cur->dt;

    return result;
}

