#include "finite-element1d.h"
#include "matrix.h"
#include "stepsize.h"
#include <math.h>
#include <stdlib.h>

double CorrectorError(struct fe1d *p, matrix *guess, matrix *soln)
{
    matrix *diff;
    int i;
    double dt, dtm1,
           d = 0, dtemp,
           denom;
    solution *s;


    s = FetchSolution(p, p->t-1);
    dt = s->dt;
    s = FetchSolution(p, p->t-2);
    dtm1 = s->dt;

    diff = mtxsub(soln, s->val);

    denom = 3*(1+dtm1/dt);

    for(i=0; i<nRows(diff); i++) {
        dtemp = fabs(val(diff, i, 0)/denom);
        if(dtemp > d)
            d = dtemp;
    }

    DestroyMatrix(diff);

    return d;
}

double StepSizeBase(struct fe1d *p, matrix *guess, matrix *soln)
{
    solution *s;
    double dt;

    s = FetchSolution(p, p->t-1);
    dt = s->dt;

    return dt * pow(ERR/fabs(CorrectorError(p, guess, soln)), 1./3.);
}

double StepSize(struct fe1d *p, matrix *guess, matrix *soln)
{
    solution *s;
    double dt_max = .0001,
           dt_cur,
           dt_new = StepSizeBase(p, guess, soln);

    s = FetchSolution(p, p->t-1);
    dt_cur = s->dt;

    if(dt_new > dt_max)
        dt_new = dt_max;

    if(dt_new > 1.3*dt_cur)
        dt_new = 1.3*dt_cur;
    //if(dt_new < .7*dt_cur)
    //    dt_new = 0.7*dt_cur;

    return dt_new;
}

double CurrentTime(struct fe1d *p, int step)
{
    int i;
    double t = 0;
    solution *s;

    for(i=0; i<step; i++) {
        s = FetchSolution(p, i);
        t += s->dt;
    }

    return t;
}
