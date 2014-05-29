/**
 * @file evap.c
 * Taylor series approximation of the rate of evaporation.
 */

/**
 * Derivative of rate of evaporation with respect to dimensionless water
 * concentration.
 * @param p Finite element problem structure
 * @param guess Current values for the dependent variables
 * @param elem Element in the mesh to calculate values for
 * @param x x-coordinate to evaluate the equation at
 * @param f1 First interpolation function to use
 * @param f2 Second interpolation function
 * @returns Value of derivative
 */
double DIDcw(struct fe1d *p, matrix *guess, Elem1D *elem, double x, int f1, int f2)
{
    double Idot = 0,
           cw0, wv0, P0, T0, phi0; /* From the current guess */

    basis *b;
    solution *s;

    b = p->b;

    /* Evaluate each of the variables at x from the guess matrix */
    s = CreateSolution(p->t, p->dt, guess);
    cw0 = EvalSoln1D(p, VARCW, elem, s, x);
    wv0 = EvalSoln1D(p, VARWV, elem, s, x);
    P0 = EvalSoln1D(p, VARP, elem, s, x);
    T0 = Tc;
    phi0 = phic;

    /* Equation */
    Idot = (evap(cw0+h, wv0, phi0, T0, P0)-evap(cw0-h, wv0, phi0, T0, P0))/(2*h);
    Idot *= b->dphi[f1](x) * b->dphi[f2](x);

    /* Clean up */
    free(s);

    return Idot;
}

/**
 * Derivative of rate of evaporation with respect to vapor mass fraction
 * @param p Finite element problem structure
 * @param guess Current values for the dependent variables
 * @param elem Element in the mesh to calculate values for
 * @param x x-coordinate to evaluate the equation at
 * @param f1 First interpolation function to use
 * @param f2 Second interpolation function
 * @returns Value of derivative
 */
double DIDwv(struct fe1d *p,
        matrix *guess,
        Elem1D *elem,
        double x, 
        int f1, int f2)
{
    double Idot = 0,
           cw0, wv0, P0, T0, phi0;

    basis *b;
    solution *s;

    b = p->b;

    /* Evaluate each of the variables at x from the guess matrix */
    s = CreateSolution(p->t, p->dt, guess);
    cw0 = EvalSoln1D(p, VARCW, elem, s, x);
    wv0 = EvalSoln1D(p, VARWV, elem, s, x);
    P0 = EvalSoln1D(p, VARP, elem, s, x);
    T0 = Tc;
    phi0 = phic;

    /* Equation */
    Idot = (evap(cw0, wv0+h, phi0, T0, P0)-evap(cw0, wv0-h, phi0, T0, P0))/(2*h);
    Idot *= b->dphi[f1](x) * b->dphi[f2](x);

    /* Clean up */
    free(s);

    return Idot;
}

/**
 * Derivative of rate of evaporation with respect to dimensionless pressure
 * @param p Finite element problem structure
 * @param guess Current values for the dependent variables
 * @param elem Element in the mesh to calculate values for
 * @param x x-coordinate to evaluate the equation at
 * @param f1 First interpolation function to use
 * @param f2 Second interpolation function
 * @returns Value of derivative
 */
double DIDP(struct fe1d *p,
        matrix *guess,
        Elem1D *elem,
        double x,
        int f1, int f2)
{
    double Idot = 0,
           cw0, wv0, P0, T0, phi0;

    basis *b;
    solution *s;

    b = p->b;

    /* Evaluate each of the variables at x from the guess matrix */
    s = CreateSolution(p->t, p->dt, guess);
    cw0 = EvalSoln1D(p, VARCW, elem, s, x);
    wv0 = EvalSoln1D(p, VARWV, elem, s, x);
    P0 = EvalSoln1D(p, VARP, elem, s, x);
    T0 = Tc;
    phi0 = phic;

    /* Equation */
    Idot = (evap(cw0, wv0, phi0, T0, P0+h)-evap(cw0, wv0, phi0, T0, P0-h))/(2*h);
    Idot *= b->dphi[f1](x) * b->dphi[f2](x);

    /* Clean up */
    free(s);

    return Idot;
}

/**
 * The portion of the Taylor series expansion that isn't a function of any of
 * the dependent variables. This is to be added to the load vector in the final
 * formulation of the problem.
 * @param p Finite element problem structure
 * @param guess Current values for the dependent variables
 * @param elem Element in the mesh to calculate values for
 * @param x x-coordinate to evaluate the equation at
 * @param f1 First interpolation function to use
 * @param f2 Second interpolation function
 * @returns Value of derivative
 */
double IdotLoad(struct fe1d *p,
        matrix *guess,
        Elem1D *elem,
        double x,
        int f1)
{
    double Idot = 0, 
           cw0, wv0, P0, T0, phi0;
    solution *s;

    /* Evaluate each of the variables at x from the guess matrix */
    s = CreateSolution(p->t, p->dt, guess);
    cw0 = EvalSoln1D(p, VARCW, elem, s, x);
    wv0 = EvalSoln1D(p, VARWV, elem, s, x);
    P0 = EvalSoln1D(p, VARP, elem, s, x);
    T0 = Tc;
    phi0 = phic;

    /* Actual equation */
    Idot = evap(cw0, wv0, phi0, T0, P0);
    Idot -= (evap(cw0+h, wv0, phi0, T0, P0)-evap(cw0-h, wv0, phi0, T0, P0))/(2*h)
        * cw0;
    Idot -= (evap(cw0, wv0+h, phi0, T0, P0)-evap(cw0, wv0+h, phi0, T0, P0))/(2*h)
        * wv0;
    Idot -= (evap(cw0, wv0, phi0, T0, P0+h)-evap(cw0, wv0, phi0, T0, P0-h))/(2*h)
        * P0;

    Idot *= p->b->phi[f1](x);

    /* Clean up */
    free(s);

    return Idot;
}

