/**
 * @file gas-phase.c
 */

/**
 * Derivative of the gas-phase residual with respect to water concentration.
 * This is used to construct the Jacobian matrix in the final formulation.
 * @param p Finite element problem structure
 * @param guess Current values for the dependent variables
 * @param elem Element in the mesh to calculate values for
 * @param x x-coordinate to evaluate the equation at
 * @param f1 First interpolation function to use
 * @param f2 Second interpolation function
 * @returns Value of derivative
 */
double DRgDcw(struct fe1d *p, matrix *guess, Elem1D *elem, double x, int f1, int f2)
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
    DKgDcw = (K_GAS(cw+h, phi, T) - K_GAS(cw-h, phi, T)) / (2*h);
    rhog = RHO_GAS(wv, T, P);
    psii = p->b->phi[f1](x);
    psij = p->b->phi[f2](x);
    dpsii = p->b->dphi[f1](x);
    dpsij = p->b->dphi[f2](x);

    /* The actual equation */
    value += -Da * rhog * DKgDcw * dpsii * psij;
    value += DIDcw(p, guess, elem, x, f1, f2);
    /* Isoparametric mapping */
    value *= IMap1D(p, elem, x);

    /* Since we don't want to get rid of our guess matrix, just delete the
     * struct we created, and not the matrix stored inside. */
    free(s);

    return value;
}

/**
 * Derivative of the gas phase residual with respect to mass fraction of vapor.
 * @param p Finite element problem structure
 * @param guess Current values for the dependent variables
 * @param elem Element in the mesh to calculate values for
 * @param x x-coordinate to evaluate the equation at
 * @param f1 First interpolation function to use
 * @param f2 Second interpolation function
 * @returns Value of derivative
 */
double DRgDwv(struct fe1d *p, matrix *guess, Elem1D *elem, double x, int f1, intf2)
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
        wv = InitValue(VARWV);
        P = InitValue(VARP);
    } else {
        cw = EvalSoln1D(p, VARCW, elem, s, x);
        wv = EvalSoln1D(p, VARWV, elem, s, x);
        P = EvalSoln1D(p, VARP, elem, s, x);
    }

    /* Calculate derivatives of Kw and D with respect to cw */
    /* TODO: Make this work with dimensionless numbers */
    DDensDwv = (RHO_GAS(wv+h, T, P) - RHO_GAS(wv-h, T, P)) / 2*h;
    psii = p->b->phi[f1](x);
    psij = p->b->phi[f2](x);
    dpsii = p->b->dphi[f1](x);
    dpsij = p->b->dphi[f2](x);

    /* The actual equation */
    value += -Da * Kg(cw, phi, T) * DDensDwv * P;
    value *= pow(dpsii, 2) * dpsij;
    value += DIDwv(p, guess, elem, x, f1, f2);
    /* Isoparametric mapping */
    value *= IMap1D(p, elem, x);

    /* Since we don't want to get rid of our guess matrix, just delete the
     * struct we created, and not the matrix stored inside. */
    free(s);

    return value;
}

/**
 * Derivative of the gas phase residual with respect to pressure.
 * @param p Finite element problem structure
 * @param guess Current values for the dependent variables
 * @param elem Element in the mesh to calculate values for
 * @param x x-coordinate to evaluate the equation at
 * @param f1 First interpolation function to use
 * @param f2 Second interpolation function
 * @returns Value of derivative
 */
double DRgDP(struct fe1d *p, matrix *guess, Elem1D *elem, double x, int f1, int f2)
{
    double value = 0, /* Final value */
           cw, /* Water concentration [-] */
           phi, /* Porosity [-] */
           T, /* Temperature [K] */
           P, /* Pressure [-] */
           Kw, /* Permeability of water [-] */
           Da, /* Darcy number [-] */
           D; /* Diffusivity [-] */
    solution *s; /* Used to get values from the current guess */

    /* Create a solution struct to store the current guessed values */
    s = CreateSolution(p->t, p->dt, guess);

    if(p->t == 0) {
        cw = InitValue(VARCW);
        wv = InitValue(VARWV);
        P = InitValue(VARP);
    } else {
        cw = EvalSoln1D(p, VARCW, elem, s, x);
        wv = EvalSoln1D(p, VARWV, elem, s, x);
        P = EvalSoln1D(p, VARP, elem, s, x);
    }

    /* Calculate derivatives of Kw and D with respect to cw */
    /* TODO: Make this work with dimensionless numbers */
    Kg = K_GAS(cw, phi, T);
    rhog = RHO_GAS(wv, T, P);
    DDensDwv = (RHO_GAS(wv+h, T, P) - RHO_GAS(wv-h, T, P)) / 2*h;
    DDensDP = (RHO_GAS(wv, T, P+h) - RHO_GAS(wv, T, P-h)) / 2*h;
    psii = p->b->phi[f1](x);
    psij = p->b->phi[f2](x);
    dpsii = p->b->dphi[f1](x);
    dpsij = p->b->dphi[f2](x);

    /* The actual equation */
    value += -Da * Kg * DDensDwv * wv * pow(dpsii, 2) * psij; //Check
    value += -Da * Kg * DDensDP * P * pow(dpsii, 2) * psij; //Check
    value += Da * rhog * Kg * dpsii * dpsij;
    value += DIDP(p, guess, elem, x, f1, f2);
    value *= IMap1D(p, elem, x);

    free(s);

    return value;
}

/**
 * Derivative of the gas phase residual (time-dependent portion) with respect to
 * water concentration. This is used in constructing the dJ matrix in the final
 * formulation.
 * @param p Finite element problem structure
 * @param guess Current values for the dependent variables
 * @param elem Element in the mesh to calculate values for
 * @param x x-coordinate to evaluate the equation at
 * @param f1 First interpolation function to use
 * @param f2 Second interpolation function
 * @returns Value of derivative
 */
double DRgDcwT(struct fe1d *p, matrix *guess, Elem1D *elem, double x, int f1, intf2)
{
    double cw,
           wv,
           phi,
           T,
           P,
           DcgDcw;
    solution *s; /* Used to get values from the current guess */
    s = CreateSolution(p->t, p->dt, guess);
    if(p->t == 0) {
        cw = InitValue(VARCW);
        wv = InitValue(VARWV);
        P = InitValue(VARP);
    } else {
        cw = EvalSoln1D(p, VARCW, elem, s, x);
        wv = EvalSoln1D(p, VARWV, elem, s, x);
        P = EvalSoln1D(p, VARP, elem, s, x);
    }
    DcgDcw = (C_GAS(cw+h, wv, phi, T, P)-C_GAS(cw-h, wv, phi, T, P))/(2*h);

    free(s);

    return DcgDcw * p->b->phi[f1](x) * p->b->phi[f2](x) * IMap1D(p, elem, x);
}

/**
 * Derivative of the gas phase residual (time-dependent portion) with respect
 * to vapor mass fraction.
 * @param p Finite element problem structure
 * @param guess Current values for the dependent variables
 * @param elem Element in the mesh to calculate values for
 * @param x x-coordinate to evaluate the equation at
 * @param f1 First interpolation function to use
 * @param f2 Second interpolation function
 * @returns Value of derivative
 */
double DRgDwvT(struct fe1d *p, matrix *guess, Elem1D *elem, double x, int f1, intf2)
{
    double cw,
           wv,
           phi,
           T,
           P,
           DcgDwv;
    solution *s; /* Used to get values from the current guess */
    s = CreateSolution(p->t, p->dt, guess);
    if(p->t == 0) {
        cw = InitValue(VARCW);
        wv = InitValue(VARWV);
        P = InitValue(VARP);
    } else {
        cw = EvalSoln1D(p, VARCW, elem, s, x);
        wv = EvalSoln1D(p, VARWV, elem, s, x);
        P = EvalSoln1D(p, VARP, elem, s, x);
    }
    DcgDwv = (C_GAS(cw, wv+h, phi, T, P)-C_GAS(cw, wv+h, phi, T, P))/(2*h);

    free(s);

    return DcgDwv * p->b->phi[f1](x) * p->b->phi[f2](x) * IMap1D(p, elem, x);
}

/**
 * Derivative of the gas phase residual (time-dependent portion) with respect
 * to pressure.
 * @param p Finite element problem structure
 * @param guess Current values for the dependent variables
 * @param elem Element in the mesh to calculate values for
 * @param x x-coordinate to evaluate the equation at
 * @param f1 First interpolation function to use
 * @param f2 Second interpolation function
 * @returns Value of derivative
 */
double DRgDPT(struct fe1d *p, matrix *guess, Elem1D *elem, double x, int f1, intf2)
{
    double cw,
           wv,
           phi,
           T,
           P,
           DcgDP;
    solution *s; /* Used to get values from the current guess */
    s = CreateSolution(p->t, p->dt, guess);
    if(p->t == 0) {
        cw = InitValue(VARCW);
        wv = InitValue(VARWV);
        P = InitValue(VARP);
    } else {
        cw = EvalSoln1D(p, VARCW, elem, s, x);
        wv = EvalSoln1D(p, VARWV, elem, s, x);
        P = EvalSoln1D(p, VARP, elem, s, x);
    }
    DcgDP = (C_GAS(cw, wv, phi, T, P+h)-C_GAS(cw, wv, phi, T, P))/(2*h);

    free(s);

    return DcgDP * p->b->phi[f1](x) * p->b->phi[f2](x) * IMap1D(p, elem, x);
}

