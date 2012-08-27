#include <stdlib.h>

#include "auxsoln.h"
#include "finite-element1d.h"
#include "matrix.h"

void FE1DInitAuxSolns(struct fe1d *p, int nvars)
{
    int i;
    p->extravars = nvars;
    p->auxsolns = (solution***) calloc(nvars, sizeof(solution**));
    for(i=0; i<nvars; i++)
        p->auxsolns[i] = (solution**) calloc(p->maxsteps, sizeof(solution*));
    return;
}

void SolveODE(struct fe1d *p, int var, int extravar,
              double (*eqn)(double, double, double), double IC)
{
    int t, i;
    int tmax = p->maxsteps;
    solution *s, *sprev, *tmp, *tmpprev;

    double dt, T, cprev;
    int rows;

    s = FetchSolution(p, 0);
    rows = mtxlen2(s->val); /* Determine how many nodes there are */

    /* Create the initial condition */
    tmp = CreateSolution(0, 0, CreateMatrix(rows, 1));
    for(i=0; i<rows; i++)
        setval(tmp->val, IC, i, 0);
    p->auxsolns[extravar][0] = tmp;

    for(t=1; t<tmax; t++) {
        sprev = s;
        s = FetchSolution(p, t);
        dt = s->dt;

        tmpprev = tmp;
        tmp = CreateSolution(t, dt, CreateMatrix(rows, 1));
        for(i=0; i<rows; i++) {
            cprev = val(tmpprev->val, i, 0);
            T = val(s->val, i, 0);

            setval(tmp->val, eqn(cprev, T, dt), i, 0);
        }
        p->auxsolns[extravar][t] = tmp;
    }

    return;
}

void PrintAuxSoln(struct fe1d *p, int var, int t)
{
    solution *s;
    s = p->auxsolns[var][t];
    mtxprnt(s->val);
}

