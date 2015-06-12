/**
 * @file predict1d.c
 * Algorithms to predict the solution at the next time step based on the
 * solution(s) at previously calculated time steps.
 */

#include <math.h>
#include <stdio.h>
#include "matrix.h"
#include "finite-element1d.h"
#include "solve.h"

/**
 * Predict the solution at the next time step by assuming that it's equal to the
 * solution at the previous time step.
 * @param p Finite element structure
 * @returns Predicted solution at t+dt
 */
matrix* PredictSolnO0(struct fe1d *p)
{
    solution *s;
    s = FetchSolution(p, p->t-1);

    /* Since the guess matrix gets deleted in the nonlinear solver, return a
     * copy of the solution. */
    return CopyMatrix(s->val);
}

/**
 * Predict the solution at the next time step using a simple first order
 * formula. This requires the solution at the previous time step and the time
 * derivative of the solution at the previous time step.
 * \f[
 * \underline{u}_{n+1}^p = \underline{u}_n + \Delta t \underline{\dot{u}}_n
 * \f]
 * @param p Finite element structure
 * @returns Predicted solution at t+dt
 */
matrix* PredictSolnO1(struct fe1d *p)
{
    solution *s; /* Solution at the previous time step */
    matrix *u, *du; /* Values for the time derivative of the solution and
                     * for the solution itself. */
    matrix *tmp1, *tmp2; /* Temporary matricies */
    double dt; /* Time step size used to calculate the previous solution. */

    /* Get the previous solution */
    s = FetchSolution(p, p->t-1);

    u = s->val;
    du = s->dval;
    dt = s->dt;

    /* If we can't use this method because the time derivative of the previous
     * solution is missing, use a less accurate one. */
    if(!du)
        return PredictSolnO0(p);

    /* u + du*dt */
    tmp1 = mtxmulconst(du, dt);
    tmp2 = mtxadd(u, tmp1);
    DestroyMatrix(tmp1); /* Delete the unneeded matrix */

    return tmp2;
}

/**
 * Predict the solution at the next time step using the second-order accurate
 * Adams-Bashforth formula. This method requires the solution at the previous
 * time step, and the time derivatives at the previous two time steps.
 * \f[
 * \underline{u}_{n+1}^p = \underline{u}_n
 * + \frac{\Delta t_n}{2}\left[\left(2+\frac{\Delta t_n}{\Delta t_{n-1}}\right)
 *      \underline{\dot{u}}_n
 * - \frac{\Delta t_n}{\Delta t_{n-1}}\underline{\dot{u}}_{n-1}\right]
 * \f]
 * @param p Finite element structure
 * @returns Predicted solution at t+dt
 */
matrix* PredictSolnO2(struct fe1d *p)
{
    /* Here, n is the last time step with a calculated solution */
    solution *sn, *sn_1; /* Solutions at t=n and t=n-1 */
    matrix *un, *dun, *dun_1; /* u(t=n), DuDt(t=n), and DuDt(t=n-1) */
    matrix *tmp1, *tmp2, *tmp3, *tmp4; /* Temporary matricies */
    double dtn, dtn_1; /* time steps at t=n-1 and t=n-2 */

    /* Get the solutions at the previous two time steps. */
    sn = FetchSolution(p, p->t-1);
    sn_1 = FetchSolution(p, p->t-2);

    /* Check to make sure we can use this method. If not, go with the first
     * order method. */
    if(!sn_1)
        return PredictSolnO1(p); /* Use a first order method if we can't */

    un = sn->val;
    dun = sn->dval;
    dun_1 = sn_1->dval;

    /* If the time derivatives aren't available, use a less accurate method. */
    if(!dun_1)
        return PredictSolnO1(p);
    if(!dun)
        return PredictSolnO0(p);

    dtn = sn->dt;
    dtn_1 = sn_1->dt;

    /* (dtn/2) * (2+dtn/dtn_1)*dun */
    tmp1 = mtxmulconst(dun, (dtn/2)*(2+dtn/dtn_1));
    /* -1 *(dtn/2) * (dtn/dtn_1)*dun_1 */
    tmp2 = mtxmulconst(dun_1, (dtn/2)*(dtn/dtn_1));
    mtxneg(tmp2);

    /* Add all three terms together */
    tmp3 = mtxadd(tmp1, tmp2);
    tmp4 = mtxadd(un, tmp3);

    /* Delete the extra matricies */
    DestroyMatrix(tmp1);
    DestroyMatrix(tmp2);
    DestroyMatrix(tmp3);

    /* Return the predicted solution */
    return tmp4;
}

