/**
 * @file vapor-phase.c
 */

/**
 * Derivative of the vapor-phase residual with respect to water concentration.
 * This is used to construct the Jacobian matrix in the final formulation.
 * @param p Finite element problem structure
 * @param guess Current values for the dependent variables
 * @param elem Element in the mesh to calculate values for
 * @param x x-coordinate to evaluate the equation at
 * @param f1 First interpolation function to use
 * @param f2 Second interpolation function
 * @returns Value of derivative
 */
double DRvDcw(struct fe1d *p, matrix *guess, Elem1D *elem, double x, int f1, int f2)
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

    /* Set the values of water concentration, mass fraction vapor, and pressure
     * at the current point and time. */
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
    DKgDcw = (K_GAS(cw+h, phi, T) - K_GAS(cw-h, phi, T)) / (2*h);
    DDDcw = (D_GAS(cw+h, T) - D_GAS(cw-h, T)) / 2*h;
    DSgDcw = (S_GAS(cw+h, phi, T) - S_GAS(cw-h, phi, T)) / (2*h);
    DxvDwv = (X_VAP(wv+h) - X_VAP(wv-h)) / (2*h);
    rhog = RHO_GAS(wv, T, P);
    Sg = S_GAS(cw, phi, T);
    Cg = CG_STAR(P);
    D = D_GAS(cw, T);
    psii = b->phi[f1](x);
    psij = b->phi[f2](x);
    dpsii = b->dphi[f1](x);
    dpsij = b->dphi[f2](x);

    /* The actual equation */
    value += -Da*rhog*DKgDcw*wv*dpsii*psii*psij;
    value += -Da*phi*Cg*Cg/rhog*DSgDcw*DxvDwv*wv*pow(dpsii, 2)*psij;
    value += -phi*Sg*Cg*Cg/rhog*DDDcw*DxvDwv*wv*pow(dpsii, 2)*psij;
    value += DIDcw(p, guess, elem, x, f1, f2);
    /* Isoparametric mapping */
    value *= IMap1D(p, elem, x);

    /* Since we don't want to get rid of our guess matrix, just delete the
     * struct we created, and not the matrix stored inside. */
    free(s);

    return value;
}

/**
 * Derivative of the vapor phase residual with respect to mass fraction of
 * vapor.
 * @param p Finite element problem structure
 * @param guess Current values for the dependent variables
 * @param elem Element in the mesh to calculate values for
 * @param x x-coordinate to evaluate the equation at
 * @param f1 First interpolation function to use
 * @param f2 Second interpolation function
 * @returns Value of derivative
 */
double DRvDwv(struct fe1d *p, matrix *guess, Elem1D *elem, double x, int f1, intf2)
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
    DrhogDwv = (RHO_GAS(wv+h, T, P) - RHO_GAS(wv-h, T, P)) / (2*h);
    DrhogDP = (RHO_GAS(wv, T, P+h) - RHO_GAS(wv, T, P-h)) / (2*h);
    DKgDcw = (K_GAS(cw+h, phi, T) - K_GAS(cw-h, phi, T)) / (2*h);
    DDDcw = (D_GAS(cw+h, T) - D_GAS(cw-h, T)) / (2*h)
    DSgDcw = (S_GAS(cw+h, phi, T) - S_GAS(cw-h, phi, T)) / (2*h);
    DxvDwv = (X_VAP(wv+h)-X_VAP(wv-h)) / (2*h);
    rhog = RHO_GAS(wv, T, P);
    Sg = S_GAS(cw, phi, T);
    Cg = CG_STAR(P);
    D = D_GAS(cw+h, T);
    psii = b->phi[f1](x);
    psij = b->phi[f2](x);
    dpsii = b->dphi[f1](x);
    dpsij = b->dphi[f2](x);

    /* The actual equation */
    value += -Da*rhog*DKgDcw*cw*dpsii*psii*psij;
    value += Da*rhog*P*pow(dpsii, 2)*psij;
    value += Da*rhog*P*psii*dpsii*dpsij;

    value += -2*Da*Kg*DrhogDwv*P*wv*pow(dpsii, 2)*psii*psij;
    value += -Da*Kg*DrhogDP*P*P*dpsii*psii*psij; // Check this line

    value += -Da*Kg*rhog*P*dpsii*dpsii*dpsij;
    value += D*phi*Cg*Cg/rhog*DSgDcw*DxvDwv*cw*dphii*dphii*dphij; // Check

    value += -phi*Sg*Cg*Cg/rhog*DDDcw*DxvDwv*cw*dpsii*dpsii*dpsij;
    value += -2*phi*Sg*D*Cg*Cg/rhog/rhog*DrhogDP*DxvDwv*P*dpsii*dpsii*dpsij;

    value += -2*phi*Sg*D*Cg*Cg/rhogD2xvDwv2*wv*dpsii*psij; // Check
    value += -phi*Sg*D*Cg/rhog*DxvDwv*dpsii*dpsij;

    value += DIDwv(p, guess, elem, x, f1, f2);
    /* Isoparametric mapping */
    value *= IMap1D(p, elem, x);

    /* Since we don't want to get rid of our guess matrix, just delete the
     * struct we created, and not the matrix stored inside. */
    free(s);

    return value;
}

/**
 * Derivative of the vapor phase residual with respect to pressure.
 * @param p Finite element problem structure
 * @param guess Current values for the dependent variables
 * @param elem Element in the mesh to calculate values for
 * @param x x-coordinate to evaluate the equation at
 * @param f1 First interpolation function to use
 * @param f2 Second interpolation function
 * @returns Value of derivative
 */
double DRvDP(struct fe1d *p, matrix *guess, Elem1D *elem, double x, int f1, int f2)
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
    DrhogDwv = (RHO_GAS(wv+h, T, P) - RHO_GAS(wv-h, T, P)) / (2*h);
    DrhogDP = (RHO_GAS(wv, T, P+h) - RHO_GAS(wv, T, P-h)) / (2*h);
    DCgDP = (CG_STAR(P+h) - CG_STAR(P-h)) / (2*h);

    /* The actual equation */
    value += Da*rhog*wv*dpsii*dpsii*psij;
    value += Da*rhog*wv*psii*dpsii*dpsij;
    value += -Da*Kg*DrhogDwv*wv*wv*dpsii*dpsii*psii*psij;

    value += -2*Da*Kg*DrhogDP*P*wv*dpsii*dpsii*psii*psij;
    value += -Da*Kg*rhog*wv*dpsii*dpsii*psij;

    value += -2*phi*Sg*D*Cg/rhog*DCgDP*DxvDwv*wv*dpsii*dpsii*dpsij;
    value += phi*Sg*D*Cg*Cg/rhog/rhog*DrhogDP*DxvDwv*wv*dpsii*dpsii*dpsij;

    value += DIDP(p, guess, elem, x, f1, f2);
    value *= IMap1D(p, elem, x);

    free(s);

    return value;
}

/**
 * Derivative of the vapor phase residual (time-dependent portion) with respect
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
double DRvDcwT(struct fe1d *p, matrix *guess, Elem1D *elem, double x, int f1, int f2)
{
    double cw,
           wv,
           phi,
           T,
           P,
           DcvDcw;
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
    DcvDcw = (C_VAP(cw+h, wv, phi, T, P)-C_VAP(cw-h, wv, phi, T, P))/(2*h);

    free(s);

    return DcvDcw * p->b->phi[f1](x) * p->b->phi[f2](x) * IMap1D(p, elem, x);
}

/**
 * Derivative of the vapor phase residual (time-dependent portion) with respect
 * to vapor mass fraction.
 * @param p Finite element problem structure
 * @param guess Current values for the dependent variables
 * @param elem Element in the mesh to calculate values for
 * @param x x-coordinate to evaluate the equation at
 * @param f1 First interpolation function to use
 * @param f2 Second interpolation function
 * @returns Value of derivative
 */
double DRvDwvT(struct fe1d *p, matrix *guess, Elem1D *elem, double x, int f1, int f2)
{
    double cw,
           wv,
           phi,
           T,
           P,
           DcvDwv;
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
    DcvDwv = (C_VAP(cw, wv+h, phi, T, P)-C_VAP(cw, wv+h, phi, T, P))/(2*h);

    free(s);

    return DcvDwv * p->b->phi[f1](x) * p->b->phi[f2](x) * IMap1D(p, elem, x);
}

/**
 * Derivative of the vapor phase residual (time-dependent portion) with respect
 * to pressure.
 * @param p Finite element problem structure
 * @param guess Current values for the dependent variables
 * @param elem Element in the mesh to calculate values for
 * @param x x-coordinate to evaluate the equation at
 * @param f1 First interpolation function to use
 * @param f2 Second interpolation function
 * @returns Value of derivative
 */
double DRvDPT(struct fe1d *p, matrix *guess, Elem1D *elem, double x, int f1, int f2)
{
    double cw,
           wv,
           phi,
           T,
           P,
           DcvDP;
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
    DcvDP = (C_VAP(cw, wv, phi, T, P+h)-C_VAP(cw, wv, phi, T, P))/(2*h);

    free(s);

    return DcgDP * p->b->phi[f1](x) * p->b->phi[f2](x) * IMap1D(p, elem, x);
}

