#include "scaling_ht.h"

scaling_ht SetupScaling(double alpha, double Tc, double Lc, double k, double h)
{
    scaling_ht stuff;

    stuff.alpha = alpha;
    stuff.Tc = Tc;
    stuff.Lc = Lc;
    stuff.k = k;
    stuff.h = h;

    return stuff;
}

double scaleTemp(scaling_ht stuff, double T)
{
    return T/stuff.Tc;
}

double uscaleTemp(scaling_ht stuff, double theta)
{
    return theta*stuff.Tc;
}

double scaleLength(scaling_ht stuff, double x)
{
    return x/stuff.Lc;
}

double uscaleLength(scaling_ht stuff, double l)
{
    return l*stuff.Lc;
}

double scaleTime(scaling_ht stuff, double t)
{
    return stuff.alpha*t/(stuff.Lc*stuff.Lc);
}

double uscaleTime(scaling_ht stuff, double XFo)
{
    return XFo*stuff.Lc*stuff.Lc/stuff.alpha;
}

double BiotNumber(scaling_ht stuff)
{
    return stuff.h*stuff.Lc/stuff.k;
}

