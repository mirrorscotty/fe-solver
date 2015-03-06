#include "finite-element1d.h"
#include "common.h"
#include "deformation.h"
#include "material-data.h"

#define CREEP CreepGina

extern choi_okos *comp_global;

double CreepGina(double t, double T, double X, double P, int deriv)
{
    double J;
    burgerse *b;
    b = CreateBurgersE();
    if(deriv)
        J = DBurgersECreep(b, t, T, X, P);
    else
        J = BurgersECreep(b, t, T, X, P);
    DestroyBurgersE(b);
    return J;
}

double CreepGinaBulk(double t, double T, double X, double P, int deriv)
{
    double nu = .45, /* Poisson ratio */
           J, /* Creep compliance */
           B; /* Bulk compliance */
    burgerse *b;
    b = CreateBurgersE();
    if(deriv)
        J = DBurgersECreep(b, t, T, X, P);
    else
        J = BurgersECreep(b, t, T, X, P);
    DestroyBurgersE(b);

    B = 6*(.5-nu)*J;
    return B;
}

double CreepLaura(double t, double T, double X, double P, int deriv)
{
    double J;
    if(deriv)
        J = DMaxwellCreepConverted(t, T, X);
    else
        J = MaxwellCreepConverted(t, T, X);
    return J;
}

double CreepLaura2(double t, double T, double X, double P, int deriv)
{
    double J;
    if(deriv)
        J = DMaxwellCreepLaura(t, T, X);
    else
        J = MaxwellCreepLaura(t, T, X);
    return J;
}

double CreepZhu(double t, double T, double X, double P, int deriv)
{
    double J, Bc = 8.2e-14;
    maxwell *m;
    m = CreateMaxwellZhu();
    if(deriv)
        J = DMaxwellCreep(m, t, T, X);
    else
        J = MaxwellCreep(m, t, T, X);
    DestroyMaxwell(m);
    return J*1e0*Bc;
}

double CreepCummings(double t, double T, double X, double P, int deriv)
{
    double J;
    if(t==0)
        t = 1e-2;
    if(deriv)
        J = DLCummingsCreep(t, T, X, P);
    else
        J = LCummingsCreep(t, T, X, P);
    return J;
}

double CreepLauraL(double t, double T, double X, double P, int deriv)
{
    double J;
    if(deriv)
        J = DLLauraCreep(t, T, X, P);
    else
        J = LLauraCreep(t, T, X, P);
    return J;
}

/* Calculate the value of the deformation gradient at the specified time and
 * spatial coordinates.
 * \f[ \underline{\underline{F}}(\underline{X}, T)\f]
 * This value is calculated from the density Choi-Okos equations, and supplied
 * to the moving mesh function to recalculate the nodal values at each time
 * step.
 * @param p Finite element problem structure
 * @param X Spatial coordinate (global)
 * @param t Time step number
 * @returns Calculated value for the deformation gradient
 */
double DeformationGrad(struct fe1d *p, double X, double t)
{
    solution *s0, *sn;
    double rho0, rhon;
    double T0 = TINIT, Tn = TINIT;
    
#ifdef CVAR
    double C0, Cn;
    choi_okos *cowet0, *cowetn;
#endif
    
    s0 = FetchSolution(p, 0);
    sn = FetchSolution(p, t);

#ifdef TVAR
    Tn = uscaleTemp(p->charvals, EvalSoln1DG(p, TVAR, sn, X, 0));
    T0 = uscaleTemp(p->charvals, EvalSoln1DG(p, TVAR, s0, X, 0));
#endif
#ifdef CVAR
    Cn = uscaleTemp(p->chardiff, EvalSoln1DG(p, CVAR, sn, X, 0));
    C0 = uscaleTemp(p->chardiff, EvalSoln1DG(p, CVAR, s0, X, 0));

    cowet0 = AddDryBasis(comp_global, C0);
    cowetn = AddDryBasis(comp_global, Cn);

    rho0 = rho(cowet0, T0);
    rhon = rho(cowetn, Tn);

    DestroyChoiOkos(cowet0);
    DestroyChoiOkos(cowetn);
#else
    rhon = rho(comp_global, Tn);
    rho0 = rho(comp_global, T0);
#endif

    return rho0/rhon;
}

double EffPorePress(double X, double T)
{
    double P = pore_press(X, T),
           P0 = pore_press(CINIT, T);
    return P;
}

double PrevStrain(struct fe1d *p, double X, int t)
{
    solution *s;
    double e = 0;
    if(t==0)
        return 0;
    s = FetchSolution(p, t);
    
    e = EvalSoln1DG(p, -1, s, X, 1);

    return e;
}

double DeformGradPc(struct fe1d *p, double X, double t)
{
#ifndef CVAR
    return 1;
#endif

    solution *s;
    double T = TINIT,
           Xdb,
           dt = p->dt,
           tf = uscaleTime(p->chardiff, t*dt),
           ti,
           e = 0,
           P,
           phi;
    int i;
 
    for(i=1; i<p->t; i++) {
        s = FetchSolution(p, i);
        ti = uscaleTime(p->chardiff, i*dt);
        Xdb = uscaleTemp(p->chardiff, EvalSoln1DG(p, CVAR, s, X, 0));
        P = EffPorePress(Xdb, T);
        e += CREEP(tf-ti, T, Xdb, -1*P, 1) * P * s->dt;
    }

    s = FetchSolution(p, t);
    Xdb = uscaleTemp(p->chardiff, EvalSoln1DG(p, CVAR, s, X, 0));
    P = EffPorePress(Xdb, T);
    e += -1*CREEP(tf, T, Xdb, -1*P, 0) * P;
    s = FetchSolution(p, 0);
    Xdb = uscaleTemp(p->chardiff, EvalSoln1DG(p, CVAR, s, X, 0));
    P = EffPorePress(Xdb, T);
    e += CREEP(0.01, T, Xdb, -1*P, 0) * P;

    //printf("X = %g, e = %g, ep = %g\n", X, e, PrevStrain(p, X, t));

    return 1-e*6*(.5-POISSON);
}

