/**
 * @file liquid-phase.c
 * Equation for calculating the fluid concentration in pasta during drying under
 * isothermal conditions with no shrinkage occuring. The one-dimensional
 * form of this equation is as follows:
 */

/**
 * Derivative of the liquid-phase residual with respect to water concentration.
 * This is used to construct the Jacobian matrix in the final formulation.
 * @param p Finite element problem structure
 * @param guess Current values for the dependent variables
 * @param elem Element in the mesh to calculate values for
 * @param x x-coordinate to evaluate the equation at
 * @param f1 First interpolation function to use
 * @param f2 Second interpolation function
 * @returns Value of derivative
 */
double DRwDcw(struct fe1d *p, matrix *guess, Elem1D *elem, double x, int f1, int f2)
{
    scaling_pasta *sc;
    sc = SetupScaling();
    double value = 0, /* Final value */
           cw, /* Water concentration [-] */
           phi, /* Porosity [-] */
           T, /* Temperature [K] */
           P, /* Pressure [-] */
           Kw, /* Permeability of water [-] */
           Da = DarcyNumber(sc), /* Darcy number [-] */
           D; /* Diffusivity [-] */
    basis *b; /* Set of interpolation functions used */
    solution *s; /* Used to get values from the current guess */

    b = p->b;

    /* Create a solution struct to store the current guessed values */
    s = CreateSolution(p->t, p->dt, guess);

    if(p->t == 0) {
        cw = InitValue(VARCW);
        P = InitValue(VARP);
    } else {
        cw = EvalSoln1D(p, VARCW, elem, s, x);
        P = EvalSoln1D(p, VARP, elem, s, x);
    }

    /* Calculate derivatives of Kw and D with respect to cw */
    /* TODO: Make this work with dimensionless numbers */
    DKwDcw = (K_WAT(cw+h, phi, T) - K_WAT(cw-h, phi, T)) / 2*h;
    DDDcw = (D_LIQ(cw+h, T) - D_LIQ(cw-h, T)) / 2*h;

    /* The actual equation */
    value = -Da*DKwDcw*P * pow(b->dphi[f1](x), 2) * b->phi[f2](x);
    value += DDDcw*cw * pow(b->dphi[f1](x), 2) * b->phi[f2](x);
    value += b->dphi[f1](x) * b->dphi[f2](x);
    value += DIDcw(p, guess, elem, x, f1, f2);
    /* Isoparametric mapping */
    value *= IMap1D(p, elem, x);

    /* Since we don't want to get rid of our guess matrix, just delete the
     * struct we created, and not the matrix stored inside. */
    free(s);

    return value;
}

/**
 * Derivative of the liquid phase residual with respect to mass fraction of
 * vapor.
 * @param p Finite element problem structure
 * @param guess Current values for the dependent variables
 * @param elem Element in the mesh to calculate values for
 * @param x x-coordinate to evaluate the equation at
 * @param f1 First interpolation function to use
 * @param f2 Second interpolation function
 * @returns Value of derivative
 */
double DRwDwv(struct fe1d *p, matrix *guess, Elem1D *elem, double x, int f1, intf2)
{
    return DIDwv(p, guess, elem, x, f1, f2;
}

/**
 * Derivative of the liquid phase residual with respect to pressure.
 * @param p Finite element problem structure
 * @param guess Current values for the dependent variables
 * @param elem Element in the mesh to calculate values for
 * @param x x-coordinate to evaluate the equation at
 * @param f1 First interpolation function to use
 * @param f2 Second interpolation function
 * @returns Value of derivative
 */
double DRwDP(struct fe1d *p, matrix *guess, Elem1D *elem, double x, int f1, int f2)
{
    double value = 0, /* Final value */
           cw, /* Water concentration [-] */
           phi, /* Porosity [-] */
           T, /* Temperature [K] */
           P, /* Pressure [-] */
           Kw, /* Permeability of water [-] */
           Da, /* Darcy number [-] */
           D; /* Diffusivity [-] */
    basis *b; /* Set of interpolation functions used */
    solution *s; /* Used to get values from the current guess */

    b = p->b;

    /* Create a solution struct to store the current guessed values */
    s = CreateSolution(p->t, p->dt, guess);

    if(p->t == 0) {
        cw = InitValue(VARCW);
        P = InitValue(VARP);
    } else {
        cw = EvalSoln1D(p, VARCW, elem, s, x);
        P = EvalSoln1D(p, VARP, elem, s, x);
    }

    /* Calculate derivatives of Kw and D with respect to cw */
    /* TODO: Make this work with dimensionless numbers */
    Kw = K_WAT(cw, phi, T);
    DKwDcw = (K_WAT(cw+h, phi, T) - K_WAT(cw-h, phi, T)) / 2*h;

    /* The actual equation */
    value = -Da*DKwDcw*cw * pow(b->dphi[f1](x), 2) * b->phi[f2](x);
    value += Da * Kw * b->dphi[f1](x) * b->dphi[f2](x);
    value += DIDP(p, guess, elem, x, f1, f2);

    value *= IMap1D(p, elem, x);

    return value;
}

/**
 * Derivative of the liquid phase residual (time-dependent portion) with respect
 * to water concentration. This is used in constructing the dJ matrix in the
 * final formulation.
 * @param p Finite element problem structure
 * @param guess Current values for the dependent variables
 * @param elem Element in the mesh to calculate values for
 * @param x x-coordinate to evaluate the equation at
 * @param f1 First interpolation function to use
 * @param f2 Second interpolation function
 * @returns Value of derivative
 */
double DRwDcwT(struct fe1d *p, matrix *guess, Elem1D *elem, double x, int f1, intf2)
{
    return p->b->dphi[f1](x) * p->b->dphi[f2](x) * IMap1D(p, elem, x);
}

/**
 * Derivative of the liquid phase residual (time-dependent portion) with respect
 * to vapor mass fraction.
 * @param p Finite element problem structure
 * @param guess Current values for the dependent variables
 * @param elem Element in the mesh to calculate values for
 * @param x x-coordinate to evaluate the equation at
 * @param f1 First interpolation function to use
 * @param f2 Second interpolation function
 * @returns Value of derivative
 */
double DRwDwvT(struct fe1d *p, matrix *guess, Elem1D *elem, double x, int f1, intf2)
{
    return 0;
}

/**
 * Derivative of the liquid phase residual (time-dependent portion) with respect
 * to pressure.
 * @param p Finite element problem structure
 * @param guess Current values for the dependent variables
 * @param elem Element in the mesh to calculate values for
 * @param x x-coordinate to evaluate the equation at
 * @param f1 First interpolation function to use
 * @param f2 Second interpolation function
 * @returns Value of derivative
 */
double DRwDPT(struct fe1d *p, matrix *guess, Elem1D *elem, double x, int f1, intf2)
{
    return 0;
}

